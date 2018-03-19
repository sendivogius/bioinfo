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

    def test_count_occurrences_approx(self):
        pattern = 'AAAAA'
        dna = 'AACAAGCTGATAAACATTTAAAGAG'
        k = 2
        self.assertEqual(11, utils.count_occurrences(dna, pattern, k))

    def test_count_occurrences_approx2(self):
        pattern = 'GAGG'
        dna = 'TTTAGAGCCTTCAGAGG'
        k = 2
        self.assertEqual(4, utils.count_occurrences(dna, pattern, k))

    def test_find_occurrences(self):
        dna = 'AAA'
        self.assertEqual([0, 1], utils.find_occurrences(dna, 'AA'))

    def test_find_occurrences2(self):
        dna = 'AGCGTCGACG'
        self.assertEqual([2, 5, 8], utils.find_occurrences(dna, 'CG'))

    def test_find_occurrences3(self):
        dna = 'GATATATGCATATACTT'
        self.assertEqual([1, 3, 9], utils.find_occurrences(dna, 'ATAT'))

    def test_find_occurrences_approx(self):
        pattern = 'ATTCTGGA'
        dna = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
        k = 3
        self.assertEqual([6, 7, 26, 27], utils.find_occurrences(dna, pattern, k))

    def test_frequent_kmers(self):
        dna = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        self.assertEqual(({'CATG', 'GCAT'}, 3), utils.frequent_kmers(dna, k))

    def test_frequent_kmers2(self):
        dna = 'AAGCAAAGGTGGG'
        k = 2
        self.assertEqual(({'AA', 'GG'}, 3), utils.frequent_kmers(dna, k))

    def test_frequent_kmers_approx(self):
        dna = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        d = 1
        self.assertEqual(({'GATG', 'ATGC', 'ATGT'}, 5), utils.frequent_kmers(dna, k, d, False))

    def test_frequent_kmers_approx_rev(self):
        dna = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        d = 1
        self.assertEqual(({'ATGT', 'ACAT'}, 9), utils.frequent_kmers(dna, k, d, True))

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

    def test_get_skew(self):
        dna = 'CATGGGCATCGGCCATACGCC'
        expected = [int(i) for i in ' 0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2'.split()]
        self.assertEqual(expected, utils.get_skew(dna))

    def test_get_min_skew_pos(self):
        dna = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
        expected = [11, 24]
        self.assertEqual(expected, utils.get_min_skew_posiion(dna))

    def test_hamming_distance(self):
        dna1 = 'GGGCCGTTGGT'
        dna2 = 'GGACCGTTGAC'
        expected = 3
        self.assertEqual(expected, utils.hamming(dna1, dna2))

    def test_immediate_neighbours(self):
        pattern = 'ACG'
        expected = {'CCG', 'TCG', 'GCG','AAG','ATG','AGG','ACA','ACC','ACT','ACG'}
        self.assertEqual(expected, utils._get_neighbourhs(pattern))

    def test_immediate_neighbours2(self):
        pattern = 'A'
        expected = {'A', 'C', 'G','T'}
        self.assertEqual(expected, utils._get_neighbourhs(pattern))

    def test_neighbours(self):
        pattern = 'AAA'
        d = 1
        expected = {'CAA', 'AAT', 'AAG', 'AAC', 'AAA', 'AGA', 'TAA', 'GAA', 'ACA', 'ATA'}
        self.assertEqual(expected, utils.get_neighbourhs(pattern, d))


    def test_get_all_kmers(self):
        pattern = 'GGACCGTTGAC'
        k = 5
        expected = {'GGACC', 'GACCG','ACCGT','CCGTT','CGTTG','GTTGA','TTGAC'}
        self.assertEqual(expected, utils.get_all_kmers(pattern, k))

if __name__ == '__main__':
    unittest.main()
