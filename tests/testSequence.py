import sys
sys.path.append('..')

import unittest
from src.sequence import Plasmid


class TestSequence(unittest.TestCase):
    
    def test_valid_input(self):
        instance = Plasmid()
        fp = 'files/valid_sequence_1.fasta'
        result = instance.parse_sequence(fp)
        self.assertEqual(result, 1)

    def test_invalid_input(self):
        instance = Plasmid()
        fp = 'files/invalid_sequence_1.fasta'
        with self.assertRaises(SystemExit) as cm:
            instance.parse_sequence(fp)

if __name__ == '__main__':
    unittest.main()