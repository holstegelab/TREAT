"""Run with the TREAT environment: python -B -m unittest discover -s tests."""

from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import functions_read_based as reads
import functions_assembly_based as assembly
from chromosome_names import resolveChromosomes, validateChromosomes


class ChromosomeNamesTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='treat_chromosomes_')
        self.root = Path(self.temp.name)
        self.sequence = 'ACGT' * 100

    def tearDown(self):
        self.temp.cleanup()

    def inputs(self, chrom, stem):
        # A misleading genome-version filename must not change naming.
        ref = self.root / (stem + '_hg19.fa')
        ref.write_text('>' + chrom + '\n' + self.sequence + '\n')
        pysam.faidx(str(ref))
        bam = self.root / (stem + '.bam')
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                  'SQ': [{'SN': chrom, 'LN': len(self.sequence)}]}
        with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
            for i in range(4):
                read = pysam.AlignedSegment()
                read.query_name = 'read_' + str(i)
                read.query_sequence = self.sequence
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 0
                read.mapping_quality = 60
                read.cigar = [(0, len(self.sequence))]
                out.write(read)
        pysam.index(str(bam))
        return str(ref), str(bam)

    def test_names_survive_extraction_and_reference_lookup(self):
        for chrom, supplied in [('4', '4'), ('chr4', 'chr4'),
                                ('GL000207.1', 'GL000207.1'),
                                ('chr4', '4'), ('4', 'chr4')]:
            with self.subTest(chrom=chrom, supplied=supplied):
                ref, bam = self.inputs(chrom, chrom + '_' + supplied)
                bed_path = self.root / (chrom + '_' + supplied + '.bed')
                bed_path.write_text(supplied + '\t41\t100\n')
                region = chrom + ':41-100'
                for module in [reads, assembly]:
                    bed, count, corrected = module.readBed(str(bed_path), str(self.root))
                    self.assertEqual(count, 1)
                    self.assertEqual(bed, {supplied: [['41', '100', supplied + ':41-100']]})
                    self.assertEqual(Path(corrected).read_text(), bed_path.read_text())
                    bed = resolveChromosomes(bed, [bam], ref, corrected)
                    self.assertEqual(bed, {chrom: [['41', '100', region]]})
                    self.assertEqual(Path(corrected).read_text(), chrom + '\t41\t100\n')
                    self.assertEqual(bed_path.read_text(), supplied + '\t41\t100\n')
                    validateChromosomes(bed, [bam], ref)
                extracted = reads.samtoolsExtract(corrected, bam, str(self.root), 'tmp_' + chrom + '.bam')
                with pysam.AlignmentFile(extracted, 'rb') as result:
                    self.assertEqual(sum(1 for _ in result), 4)
                ref_reads, fasta = reads.measureDistance_reference(corrected, 10, ref, str(self.root))
                self.assertEqual(ref_reads[0][1], region)
                self.assertEqual(ref_reads[0][-4], self.sequence[40:100])
                ref_asm = assembly.measureDistance_reference_opt([region], ref, 10)
                self.assertEqual(ref_asm[0][1], region)
                self.assertEqual(ref_asm[0][-4], self.sequence[40:100])
                self.assertEqual(reads.checkIntervals(bed, chrom, 0, 400, 10), [region])

    def test_mismatch_is_reported_for_fasta_and_every_bam(self):
        plain_ref, plain_bam = self.inputs('4', 'plain')
        prefixed_ref, prefixed_bam = self.inputs('chr4', 'prefixed')
        bed = {'4': [['41', '100', '4:41-100']]}
        with self.assertRaisesRegex(ValueError, 'reference FASTA.*missing BED chromosome name'):
            validateChromosomes(bed, [plain_bam], prefixed_ref)
        with self.assertRaisesRegex(ValueError, 'BAM.*prefixed.bam.*missing BED chromosome name'):
            validateChromosomes(bed, [plain_bam, prefixed_bam], plain_ref)
        with self.assertRaisesRegex(ValueError, 'must match exactly'):
            validateChromosomes({'chr4': []}, [plain_bam], plain_ref)
        corrected = self.root / 'corrected.bed'
        corrected.write_text('4\t41\t100\n')
        for ref, bams in [(prefixed_ref, [plain_bam]),
                          (plain_ref, [plain_bam, prefixed_bam])]:
            with self.assertRaisesRegex(ValueError, 'Cannot resolve BED chromosome 4'):
                resolveChromosomes(bed, bams, ref, str(corrected))
            self.assertEqual(corrected.read_text(), '4\t41\t100\n')
        with self.assertRaisesRegex(ValueError, 'Cannot resolve BED chromosome unknown'):
            resolveChromosomes({'unknown': []}, [plain_bam], plain_ref, str(corrected))

    def test_exact_names_take_priority_when_both_exist(self):
        ref = self.root / 'both.fa'
        ref.write_text('>4\n' + self.sequence + '\n>chr4\n' + self.sequence + '\n')
        pysam.faidx(str(ref))
        bam = self.root / 'both.bam'
        with pysam.AlignmentFile(str(bam), 'wb', header={
                'SQ': [{'SN': chrom, 'LN': 400} for chrom in ['4', 'chr4']]}) as out:
            pass
        bed = {'4': [['41', '100', '4:41-100']],
               'chr4': [['41', '100', 'chr4:41-100']]}
        path = self.root / 'exact.bed'
        path.write_text('4\t41\t100\nchr4\t41\t100\n')
        self.assertEqual(resolveChromosomes(bed, [str(bam)], str(ref), str(path)), bed)

    def test_output_preserves_reference_annotation(self):
        for chrom in ['4', 'chr4', 'GL000207.1']:
            with self.subTest(chrom=chrom):
                region = chrom + ':41-100'
                # Empty reads calls still need the correct reference annotation.
                columns = ['READ_NAME', 'HAPLOTAG', 'REGION', 'PASSES', 'READ_QUALITY',
                           'LEN_SEQUENCE_FOR_TRF', 'START_TRF', 'END_TRF', 'type',
                           'SAMPLE_NAME', 'POLISHED_HAPLO', 'DEPTH', 'CONSENSUS_MOTIF',
                           'CONSENSUS_MOTIF_COPIES', 'MOTIF_REF', 'REFERENCE_MOTIF_COPIES',
                           'SEQUENCE_WITH_PADDINGS', 'SEQUENCE_FOR_TRF']
                record, sequences = reads.prepareOutputs(pd.DataFrame(columns=columns),
                    {region: ['ACGT', 60, 15]}, region, 'reads', [0, 0])
                self.assertEqual(record[0], chrom)
                self.assertEqual(record[3], 60)
                self.assertEqual(record[7], 'ACGT;15')
                seq = self.sequence[40:100]
                data = pd.DataFrame([{'SAMPLE': 'sample', 'REGION': region,
                    'HAPLOTYPE': 1, 'SEQUENCE': seq, 'POLISHED_HAPLO': 60,
                    'CONSENSUS_MOTIF': 'ACGT', 'CONSENSUS_MOTIF_COPIES': 15,
                    'COVERAGE_HAPLO': 4}])
                record = assembly.prepareOutputs_opt([region], data,
                    {region: ['ACGT', 15, 60, seq]}, ['sample'])[0]
                self.assertEqual(record[0], chrom)
                self.assertEqual(record[3], seq)
                self.assertEqual(record[7], 'ACGT;15;60')
                self.assertEqual(record[-1].split(';')[1], '0|0')


if __name__ == '__main__':
    unittest.main()
