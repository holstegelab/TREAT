"""Regression tests for the flanking sequence reported in GitHub issue #6."""

from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import functions_assembly_based as assembly


class AssemblyTrimmingTest(unittest.TestCase):
    def test_issue6_sequence_matches_reference_with_anchor(self):
        # Exact GRCh38 sequence at chr10:16854109-16854118 with otter's
        # 20-base flanks. The old default window=10 retained those flanks.
        padded = 'TGAGGATCTCAACTGTTCCAAGGAGGAGGTGATTTGGATCTGGGTCTTG'
        self.assertEqual(padded[9:-10], 'CAACTGTTCCAAGGAGGAGGTGATTTGGAT')
        with tempfile.TemporaryDirectory(prefix='treat_issue6_') as temp:
            root = Path(temp)
            fasta = root / 'sample.fa'
            fasta.write_text('>sample.fa#chr10:16854109-16854118#0'
                             '#tc:i:30#ac:i:30#sc:i:30\n' + padded + '\n')
            results = assembly.run_trf_asm_opt(str(fasta), 20)
            self.assertTrue(results)
            for row in results:
                self.assertEqual(row[15], 'AAGGAGGAGG')
                self.assertEqual(row[16], 10)
            # Move the same locus to 20-29 in a tiny reference, retaining
            # the same preceding anchor base and target sequence.
            reference = root / 'reference.fa'
            reference.write_text('>chr10\n' + padded + '\n')
            ref_rows = assembly.measureDistance_reference_opt(
                ['chr10:20-29'], str(reference), 10)
            self.assertEqual(ref_rows[0][-4], results[0][15])
            self.assertEqual(padded[20:29], 'AGGAGGAGG')


if __name__ == '__main__':
    unittest.main()
