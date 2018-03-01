import unittest
from . import utils


class UtilsTests(unittest.TestCase):
    def test_reverse_strand(self):
        dna = 'AAAACCCGGT'
        reversed = 'ACCGGGTTTT'
        self.assertEqual(reversed, utils.reverse_strand(dna))

    def test_count_occurrences(self):
        dna = 'AAA'
        self.assertEqual(2, utils.count_occurrences(dna, 'AA'))

    def test_count_occurrences2(self):
        dna = 'AGCGTGACG'
        self.assertEqual(2, utils.count_occurrences(dna, 'CG'))

    def test_find_occurrences(self):
        dna = 'AAA'
        self.assertEqual([0, 1], utils.find_occurrences(dna, 'AA'))

    def test_find_occurrences2(self):
        dna = 'AGCGTCGACG'
        self.assertEqual([2, 5, 8], utils.find_occurrences(dna, 'CG'))

    def test_find_occurrences3(self):
        dna = 'GATATATGCATATACTT'
        self.assertEqual([1, 3, 9], utils.find_occurrences(dna, 'ATAT'))

    def test_frequent_kmers(self):
        dna = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        self.assertEqual(({'CATG', 'GCAT'}, 3), utils.frequent_kmers(dna, k))

    def test_frequent_kmers2(self):
        dna = 'AAGCAAAGGTGGG'
        k = 2
        self.assertEqual(({'AA', 'GG'}, 3), utils.frequent_kmers(dna, k))

    def test_find_clumps(self):
        dna = 'CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
        k = 5
        L = 50
        t = 4
        self.assertEqual({'CGACA', 'GAAGA'}, utils.find_clumps(dna, k, L, t))

    def test_pattern_to_number(self):
        dna = 'ATGCAA'
        self.assertEqual(912, utils.pattern_to_number(dna))

    def test_pattern_to_number2(self):
        dna = 'CCCATTC'
        self.assertEqual(5437, utils.pattern_to_number(dna))

    def test_pattern_to_number3(self):
        dna = 'TAGTTCCATCGCAGAG'
        self.assertEqual(3419724066, utils.pattern_to_number(dna))

    def test_number_to_pattern(self):
        number = 5437
        length = 8
        self.assertEqual('ACCCATTC', utils.number_to_pattern(number, length))

    def test_number_to_pattern2(self):
        number = 45
        length = 4
        self.assertEqual('AGTC', utils.number_to_pattern(number, length))

    def test_number_to_pattern3(self):
        number = 8111
        length = 9
        self.assertEqual('AACTTGGTT', utils.number_to_pattern(number, length))

    def test_compute_frequencies(self):
        dna = 'ACGCGGCTCTGAAA'
        k = 2
        expected = [int(i) for i in '2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0'.split()]
        self.assertEqual(expected, utils.compute_frequencies(dna, k))



if __name__ == '__main__':
    unittest.main()
