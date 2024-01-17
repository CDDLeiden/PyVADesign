import sys
sys.path.append('..')

import unittest
from src.mutation import Mutation


class TestMutation(unittest.TestCase):
    
    def test_valid_input(self):
        instance = Mutation()
        fp = 'files/valid_mutations_1.txt'
        result = instance.parse_mutations(fp)
        self.assertEqual(result, 1)

    def test_invalid_input(self):
        instance = Mutation()
        fp = 'files/invalid_mutations_1.txt'
        with self.assertRaises(SystemExit) as cm:
            instance.parse_mutations(fp)

if __name__ == '__main__':
    unittest.main()