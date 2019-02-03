import unittest

from utils import is_dna


class UtilsTests(unittest.TestCase):
    def test_is_dna_empty(self):
        self.assertTrue(is_dna(''))

    def test_is_dna_valid(self):
        self.assertTrue(is_dna('ACTG'))

    def test_is_dna_invalid(self):
        self.assertFalse(is_dna('ACTGD'))


if __name__ == '__main__':
    unittest.main()
