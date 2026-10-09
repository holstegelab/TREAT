"""Regression tests for empty assemblies and malformed/failed TRF output."""

import gzip
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import functions_assembly_based as assembly
import functions_read_based as reads


class TRFParsingTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='treat_trf_')
        self.fasta = Path(self.temp.name) / 'sample.fa'
        self.fasta.write_text('>chr4:41-100;sample;read1\n' + 'AC' * 30 + '\n')
        self.row = '1 60 2 30.0 2 100 0 120 50 50 0 0 1.00 AC ' + 'AC' * 30 + ' . .'

    def tearDown(self):
        self.temp.cleanup()

    def parse(self, stdout='', stderr='', returncode=0):
        result = subprocess.CompletedProcess([], returncode, stdout, stderr)
        with patch.object(reads.subprocess, 'run', return_value=result):
            return reads.run_trf(0, [str(self.fasta)], 'reads')

    def test_valid_results_and_multiple_headers(self):
        results = self.parse('@chr4:41-100;sample;read1\n' + self.row +
                             '\n@chr4:41-100;sample;read2\n' + self.row + '\n')
        self.assertEqual([row[0] for row in results],
                         ['read1_chr4:41-100', 'read2_chr4:41-100'])
        self.assertTrue(all(len(row) == 19 for row in results))

    def test_empty_input_and_successful_no_matches(self):
        for output in ['', '@chr4:41-100;sample;read1\n']:
            self.assertEqual(self.parse(output), [['NA'] * 19])
        self.fasta.write_text(' \n\t\n')
        with patch.object(reads.subprocess, 'run') as run:
            self.assertEqual(reads.run_trf(0, [str(self.fasta)], 'reads'), [['NA'] * 19])
            run.assert_not_called()

    def test_diagnostics_before_header_and_bad_header(self):
        with self.assertRaisesRegex(ValueError, 'before a sequence header'):
            self.parse('Error: Could not load sequence. Empty file or bad format.\n')
        with self.assertRaisesRegex(ValueError, 'Malformed TRF sequence header'):
            self.parse('@invalid\n' + self.row)

    def test_bad_result_shape_and_numeric_fields(self):
        for row in ['1 2 3', self.row.replace('1 60', 'invalid 60', 1)]:
            with self.assertRaisesRegex(ValueError, 'Malformed TRF result'):
                self.parse('@chr4:41-100;sample;read1\n' + row)

    def test_command_failure_includes_diagnostics(self):
        with self.assertRaisesRegex(RuntimeError, 'exit code 2.*Cannot read FASTA'):
            self.parse(stderr='Cannot read FASTA', returncode=2)

    @unittest.skipUnless(shutil.which('trf'), 'TRF executable required')
    def test_actual_trf(self):
        result = reads.run_trf(0, [str(self.fasta)], 'reads')
        self.assertEqual(result[0][0], 'read1_chr4:41-100')
        self.assertEqual(result[0][15], 'AC')


class EmptyAssemblyTest(unittest.TestCase):
    def test_empty_fasta_returns_no_annotations(self):
        with tempfile.TemporaryDirectory(prefix='treat_empty_assembly_') as temp:
            fasta = Path(temp) / 'empty.fa'
            for content in ['', ' \n\t\n']:
                fasta.write_text(content)
                self.assertEqual(assembly.run_trf_asm_opt(str(fasta), 20), [])

    def test_otter_failure_is_not_a_missing_assembly(self):
        with tempfile.TemporaryDirectory(prefix='treat_failed_otter_') as temp:
            (Path(temp) / 'otter_local_asm').mkdir()
            with patch.object(assembly.subprocess, 'run', return_value=
                              subprocess.CompletedProcess([], 2)):
                with self.assertRaisesRegex(RuntimeError, 'otter failed.*exit code 2'):
                    assembly.assembly_otter_opt('sample.bam', temp, 'ref.fa', 'target.bed',
                                               1, 20, 'False', 200, '500,0.1', 0.9, 2)

    @unittest.skipUnless(shutil.which('otter') and shutil.which('samtools'),
                         'otter and samtools required')
    def test_full_pipeline_keeps_empty_samples_and_unassembled_regions(self):
        with tempfile.TemporaryDirectory(prefix='treat_assembly_pipeline_') as temp:
            root = Path(temp)
            reference = root / 'reference.fa'
            sequence = 'ACGT' * 100
            reference.write_text('>chr4\n' + sequence + '\n')
            pysam.faidx(str(reference))
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                      'SQ': [{'SN': 'chr4', 'LN': 400}]}
            bams = []
            for name, coverage in [('empty', 0), ('covered', 30)]:
                bam = root / (name + '.bam')
                with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
                    for i in range(coverage):
                        read = pysam.AlignedSegment()
                        read.query_name = 'read_' + str(i)
                        read.query_sequence = sequence[:150]
                        read.flag = 0
                        read.reference_id = 0
                        read.reference_start = 0
                        read.mapping_quality = 60
                        read.cigar = [(0, 150)]
                        out.write(read)
                pysam.index(str(bam))
                bams.append(str(bam))
            bed_path = root / 'target.bed'
            bed_path.write_text('chr4\t41\t100\nchr4\t201\t260\nchr4\t301\t360\n')
            for bam_list, label in [([bams[0]], 'all_empty'), (bams, 'mixed')]:
                output = root / label
                output.mkdir()
                bed, count, corrected = assembly.readBed(str(bed_path), str(output))
                data = assembly.otterPipeline_opt(str(output), 1, str(reference), corrected,
                        bam_list, count, 20, 10, bed, 'False', 200, '500,0.1', 0.9, 2)
                assembly.haplotyping_steps_opt(data, 1, 'otter', str(output), bam_list)
                with gzip.open(str(output / 'sample.vcf.gz'), 'rt') as vcf:
                    lines = [line.rstrip().split('\t') for line in vcf
                             if not line.startswith('##')]
                self.assertEqual(lines[0][9:], [Path(bam).stem for bam in bam_list])
                self.assertEqual(len(lines[1:]), 3)
                for record in lines[1:]:
                    self.assertEqual(record[9], 'NO_ASSEMBLY;.|.;.|.;.|.;.|.;.|.;.|.')
                    if label == 'mixed':
                        self.assertEqual(record[10].split(';')[0],
                                         'PASS' if record[2] == 'chr4:41-100' else 'NO_ASSEMBLY')


if __name__ == '__main__':
    unittest.main()
