import unittest
import motifs


class MotifsTests(unittest.TestCase):
    def test_reverse_complement_strand(self):
        dna = 'AAAACCCGGT'
        reversed = 'ACCGGGTTTT'
        actual = motifs.reverse_complement_strand(dna)
        self.assertEqual(reversed, actual)

    def test_count_occurrences_overlapping(self):
        dna = 'AAA'
        expected = 2
        actual = motifs.count_occurrences(dna, 'AA')
        self.assertEqual(expected, actual)

    def test_count_occurrences_non_overlapping(self):
        dna = 'AGCGTGACG'
        expected = 2
        actual = motifs.count_occurrences(dna, 'CG')
        self.assertEqual(expected, actual)

    def test_count_occurrences_hamming_1(self):
        dna = 'AACAAGCTGATAAACATTTAAAGAG'
        pattern = 'AAAAA'
        k = 1
        expected = 4  # AACAA, ATAAA, AAACA, AAAGA
        actual = motifs.count_occurrences(dna, pattern, k)
        self.assertEqual(expected, actual)

    def test_count_occurrences_hamming_2(self):
        pattern = 'GAGG'
        dna = 'TTTAGAGCCTTCAGAGG'
        k = 2
        expected = 4  # TAGA, GAGC, CAGA, GAGG
        actual = motifs.count_occurrences(dna, pattern, k)
        self.assertEqual(expected, actual)

    def test_find_occurrences_overlapping(self):
        dna = 'AAA'
        expected = [0, 1]
        actual = motifs.find_occurrences(dna, 'AA')
        self.assertEqual(expected, actual)

    def test_find_occurrences_non_overlapping(self):
        dna = 'AGCGTCGACG'
        expected = [2, 5, 8]
        actual = motifs.find_occurrences(dna, 'CG')
        self.assertEqual(expected, actual)

    def test_find_occurrences_partially_overlapping(self):
        dna = 'GATATATGCATATACTT'
        expected = [1, 3, 9]
        actual = motifs.find_occurrences(dna, 'ATAT')
        self.assertEqual(expected, actual)

    def test_find_occurrences_hamming(self):
        dna = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
        pattern = 'ATTCTGGA'
        k = 3
        expected = [6, 7, 26, 27]
        actual = motifs.find_occurrences(dna, pattern, k)
        self.assertEqual(expected, actual)

    def test_frequent_kmers(self):
        dna = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        expected = {'CATG', 'GCAT'}, 3
        actual = motifs.frequent_kmers(dna, k)
        self.assertEqual(expected, actual)

    def test_frequent_kmers2(self):
        dna = 'AAGCAAAGGTGGG'
        k = 2
        expected = {'AA', 'GG'}, 3
        actual = motifs.frequent_kmers(dna, k)
        self.assertEqual(expected, actual)

    def test_frequent_kmers_hamming(self):
        dna = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        d = 1
        expected = {'GATG', 'ATGC', 'ATGT'}, 5
        actual = motifs.frequent_kmers(dna, k, d, False)
        self.assertEqual(expected, actual)

    def test_frequent_kmers_hamming_reverse(self):
        dna = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        d = 1
        expected = {'ATGT', 'ACAT'}, 9
        actual = motifs.frequent_kmers(dna, k, d, True)
        self.assertEqual(expected, actual)

    def test_find_clumps(self):
        dna = 'CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
        k = 5
        L = 50
        t = 4
        expected = {'CGACA', 'GAAGA'}
        actual = motifs.find_clumps(dna, k, L, t)
        self.assertEqual(expected, actual)

    def test_pattern_to_number(self):
        dna = 'ATGCAA'
        expected = 912
        actual = motifs.pattern_to_number(dna)
        self.assertEqual(expected, actual)

    def test_pattern_to_number2(self):
        dna = 'CCCATTC'
        expected = 5437
        actual = motifs.pattern_to_number(dna)
        self.assertEqual(expected, actual)

    def test_pattern_to_number3(self):
        dna = 'TAGTTCCATCGCAGAG'
        expected = 3419724066
        actual = motifs.pattern_to_number(dna)
        self.assertEqual(expected, actual)

    def test_number_to_pattern(self):
        number = 5437
        length = 8
        expected = 'ACCCATTC'
        actual = motifs.number_to_pattern(number, length)
        self.assertEqual(expected, actual)

    def test_number_to_pattern2(self):
        number = 45
        length = 4
        expected = 'AGTC'
        actual = motifs.number_to_pattern(number, length)
        self.assertEqual(expected, actual)

    def test_number_to_pattern3(self):
        number = 8111
        length = 9
        expected = 'AACTTGGTT'
        actual = motifs.number_to_pattern(number, length)
        self.assertEqual(expected, actual)

    def test_get_skew(self):
        dna = 'CATGGGCATCGGCCATACGCC'
        expected = [int(i) for i in ' 0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2'.split()]
        actual = motifs.get_skew(dna)
        self.assertEqual(expected, actual)

    def test_get_min_skew_pos(self):
        dna = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
        expected = [11, 24]
        actual = motifs.get_min_skew_position(dna)
        self.assertEqual(expected, actual)

    def test_hamming_distance(self):
        dna1 = 'GGGCCGTTGGT'
        dna2 = 'GGACCGTTGAC'
        expected = 3
        actual = motifs.hamming(dna1, dna2)
        self.assertEqual(expected, actual)

    def test_immediate_neighbours(self):
        pattern = 'ACG'
        expected = {'CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG', 'ACA', 'ACC', 'ACT', 'ACG'}
        actual = motifs._get_neighbours(pattern)
        self.assertEqual(expected, actual)

    def test_immediate_neighbours_one_nucletide(self):
        pattern = 'A'
        expected = {'A', 'C', 'G', 'T'}
        actual = motifs._get_neighbours(pattern)
        self.assertEqual(expected, actual)

    def test_neighbours_d1(self):
        pattern = 'AAA'
        d = 1
        expected = {'CAA', 'AAT', 'AAG', 'AAC', 'AAA', 'AGA', 'TAA', 'GAA', 'ACA', 'ATA'}
        actual = motifs.get_neighbours(pattern, d)
        self.assertEqual(expected, actual)

    def test_neighbours_d2(self):
        pattern = 'ACT'
        d = 2
        expected = {'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC',
                    'ATG', 'ATT', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGT', 'CTT', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
                    'GGT', 'GTT', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGT', 'TTT'}
        actual = motifs.get_neighbours(pattern, d)
        self.assertEqual(expected, actual)

    def test_get_all_kmers(self):
        pattern = 'GGACCGTTGAC'
        k = 5
        expected = {'GGACC', 'GACCG', 'ACCGT', 'CCGTT', 'CGTTG', 'GTTGA', 'TTGAC'}
        actual = motifs.get_all_kmers(pattern, k)
        self.assertEqual(expected, actual)

    def test_get_profile_no_pseudo_no_relative(self):
        patterns = ['TCGGGGGTTTTT', 'CCGGTGACTTAC', 'ACGGGGATTTTC', 'TTGGGGACTTTT', 'AAGGGGACTTCC', 'TTGGGGACTTCC',
                    'TCGGGGATTCAT', 'TCGGGGATTCCT', 'TAGGGGAACTAC', 'TCGGGTATAACC']
        expected = [{'A': 2, 'C': 1, 'G': 0, 'T': 7},
                    {'A': 2, 'C': 6, 'G': 0, 'T': 2},
                    {'A': 0, 'C': 0, 'G': 10, 'T': 0},
                    {'A': 0, 'C': 0, 'G': 10, 'T': 0},
                    {'A': 0, 'C': 0, 'G': 9, 'T': 1},
                    {'A': 0, 'C': 0, 'G': 9, 'T': 1},
                    {'A': 9, 'C': 0, 'G': 1, 'T': 0},
                    {'A': 1, 'C': 4, 'G': 0, 'T': 5},
                    {'A': 1, 'C': 1, 'G': 0, 'T': 8},
                    {'A': 1, 'C': 2, 'G': 0, 'T': 7},
                    {'A': 3, 'C': 4, 'G': 0, 'T': 3},
                    {'A': 0, 'C': 6, 'G': 0, 'T': 4}, ]
        actual = motifs.get_profile(patterns, relative=False, pseudocounts=False)
        self.assertEqual(expected, actual)

    def test_get_profile_pseudo_no_relative(self):
        patterns = ['TCGGGGGTTTTT', 'CCGGTGACTTAC', 'ACGGGGATTTTC', 'TTGGGGACTTTT', 'AAGGGGACTTCC', 'TTGGGGACTTCC',
                    'TCGGGGATTCAT', 'TCGGGGATTCCT', 'TAGGGGAACTAC', 'TCGGGTATAACC']
        expected = [{'A': 3, 'C': 2, 'G': 1, 'T': 8},
                    {'A': 3, 'C': 7, 'G': 1, 'T': 3},
                    {'A': 1, 'C': 1, 'G': 11, 'T': 1},
                    {'A': 1, 'C': 1, 'G': 11, 'T': 1},
                    {'A': 1, 'C': 1, 'G': 10, 'T': 2},
                    {'A': 1, 'C': 1, 'G': 10, 'T': 2},
                    {'A': 10, 'C': 1, 'G': 2, 'T': 1},
                    {'A': 2, 'C': 5, 'G': 1, 'T': 6},
                    {'A': 2, 'C': 2, 'G': 1, 'T': 9},
                    {'A': 2, 'C': 3, 'G': 1, 'T': 8},
                    {'A': 4, 'C': 5, 'G': 1, 'T': 4},
                    {'A': 1, 'C': 7, 'G': 1, 'T': 5}, ],
        actual = motifs.get_profile(patterns, relative=False, pseudocounts=True)
        self.assertEqual(expected, actual)

    def test_get_profile_no_pseudo_relative(self):
        patterns = ['TCGGGGGTTTTT', 'CCGGTGACTTAC', 'ACGGGGATTTTC', 'TTGGGGACTTTT', 'AAGGGGACTTCC', 'TTGGGGACTTCC',
                    'TCGGGGATTCAT', 'TCGGGGATTCCT', 'TAGGGGAACTAC', 'TCGGGTATAACC']
        expected = [{'A': .2, 'C': .1, 'G': 0, 'T': .7},
                    {'A': .2, 'C': .6, 'G': 0, 'T': .2},
                    {'A': 0, 'C': 0, 'G': 1, 'T': 0},
                    {'A': 0, 'C': 0, 'G': 1, 'T': 0},
                    {'A': 0, 'C': 0, 'G': .9, 'T': .1},
                    {'A': 0, 'C': 0, 'G': .9, 'T': .1},
                    {'A': .9, 'C': 0, 'G': .1, 'T': 0},
                    {'A': .1, 'C': .4, 'G': 0, 'T': .5},
                    {'A': .1, 'C': .1, 'G': 0, 'T': .8},
                    {'A': .1, 'C': .2, 'G': 0, 'T': .7},
                    {'A': .3, 'C': .4, 'G': 0, 'T': .3},
                    {'A': 0, 'C': .6, 'G': 0, 'T': .4}]
        actual = motifs.get_profile(patterns, relative=True, pseudocounts=False)
        self.assertEqual(expected, actual)

    def test_get_profile_pseudo_relative(self):
        patterns = ['TCGGGGGTTTTT', 'CCGGTGACTTAC', 'ACGGGGATTTTC', 'TTGGGGACTTTT', 'AAGGGGACTTCC', 'TTGGGGACTTCC',
                    'TCGGGGATTCAT', 'TCGGGGATTCCT', 'TAGGGGAACTAC', 'TCGGGTATAACC']
        expected = [{'A': 3. / 14, 'C': 2. / 14, 'G': 1. / 14, 'T': 8. / 14},
                    {'A': 3. / 14, 'C': 7. / 14, 'G': 1. / 14, 'T': 3. / 14},
                    {'A': 1. / 14, 'C': 1. / 14, 'G': 11. / 14, 'T': 1. / 14},
                    {'A': 1. / 14, 'C': 1. / 14, 'G': 11. / 14, 'T': 1. / 14},
                    {'A': 1. / 14, 'C': 1. / 14, 'G': 10. / 14, 'T': 2. / 14},
                    {'A': 1. / 14, 'C': 1. / 14, 'G': 10. / 14, 'T': 2. / 14},
                    {'A': 10. / 14, 'C': 1. / 14, 'G': 2. / 14, 'T': 1. / 14},
                    {'A': 2. / 14, 'C': 5. / 14, 'G': 1. / 14, 'T': 6. / 14},
                    {'A': 2. / 14, 'C': 2. / 14, 'G': 1. / 14, 'T': 9. / 14},
                    {'A': 2. / 14, 'C': 3. / 14, 'G': 1. / 14, 'T': 8. / 14},
                    {'A': 4. / 14, 'C': 5. / 14, 'G': 1. / 14, 'T': 4. / 14},
                    {'A': 1. / 14, 'C': 7. / 14, 'G': 1. / 14, 'T': 5. / 14}, ]
        actual = motifs.get_profile(patterns, relative=True, pseudocounts=True)
        self.assertEqual(expected, actual)

    def test_get_consensus_string(self):
        patterns = ['TCGGGGGTTTTT', 'CCGGTGACTTAC', 'ACGGGGATTTTC', 'TTGGGGACTTTT', 'AAGGGGACTTCC', 'TTGGGGACTTCC',
                    'TCGGGGATTCAT', 'TCGGGGATTCCT', 'TAGGGGAACTAC', 'TCGGGTATAACC']
        expected = 'TCGGGGATTTCC'
        actual = motifs.get_consensus_string(patterns)
        self.assertEqual(expected, actual)

    def test_score_motif_mismatches(self):
        patterns = ['TCGGGGGTTTTT', 'CCGGTGACTTAC', 'ACGGGGATTTTC', 'TTGGGGACTTTT', 'AAGGGGACTTCC', 'TTGGGGACTTCC',
                    'TCGGGGATTCAT', 'TCGGGGATTCCT', 'TAGGGGAACTAC', 'TCGGGTATAACC']
        expected = 30
        actual = motifs.score_motifs(patterns)
        self.assertEqual(expected, actual)

    def test_score_motif_entropy(self):
        patterns = ['TCGGGGGTTTTT', 'CCGGTGACTTAC', 'ACGGGGATTTTC', 'TTGGGGACTTTT', 'AAGGGGACTTCC', 'TTGGGGACTTCC',
                    'TCGGGGATTCAT', 'TCGGGGATTCCT', 'TAGGGGAACTAC', 'TCGGGTATAACC']
        expected = 9.9163
        actual = motifs.score_motifs(patterns, entropy=True)
        self.assertAlmostEqual(expected, actual, places=4)

    def test_dist_pattern_single_dna(self):
        pattern = 'AAA'
        dna = 'TTACCTTAAC'
        expected = 1  # for AAC
        actual = motifs.dist(pattern, dna)
        self.assertEqual(expected, actual)

    def test_dist_pattern_dna_list(self):
        pattern = 'AAA'
        dnas = ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT']
        expected = 5
        actual = motifs.dist(pattern, dnas)
        self.assertEqual(expected, actual)

    def test_median_string(self):
        patterns = ['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTTCGGGACAG']
        expected = ('GAC', 2)
        actual = motifs.median_string(patterns, 3)
        self.assertEqual(expected, actual)

    def test_score_kmer_on_profile(self):
        profile = [{'A': .2, 'C': .4, 'G': .3, 'T': .1},
                   {'A': .2, 'C': .3, 'G': .3, 'T': .2},
                   {'A': .3, 'C': .1, 'G': .5, 'T': .1},
                   {'A': .2, 'C': .5, 'G': .2, 'T': .1},
                   {'A': .3, 'C': .1, 'G': .4, 'T': .2}]
        kmer = 'CCGAG'
        expected = .4 * .3 * .5 * .2 * .4
        actual = motifs.score_kmer_on_profile(kmer, profile)
        self.assertAlmostEqual(expected, actual)

    def test_profile_most_probable_string(self):
        profile = [{'A': .2, 'C': .4, 'G': .3, 'T': .1},
                   {'A': .2, 'C': .3, 'G': .3, 'T': .2},
                   {'A': .3, 'C': .1, 'G': .5, 'T': .1},
                   {'A': .2, 'C': .5, 'G': .2, 'T': .1},
                   {'A': .3, 'C': .1, 'G': .4, 'T': .2}]
        dna = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
        expected = 'CCGAG'
        actual = motifs.most_probable_kmer_from_profile(dna, profile)
        self.assertEqual(expected, actual)

    def test_greedy_motif_search(self):
        dna = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
        k = 3
        expected = ['TTC', 'ATC', 'TTC', 'ATC', 'TTC'], 2
        actual = motifs.greedy_motif_search(dna, k)
        self.assertEqual(expected, actual)

    def test_randomized_motif_search(self):
        dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
               'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
               'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
               'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
               'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
        k = 8
        expected = ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TTCAGGTG', 'TCCACGTG']
        expected_score = motifs.score_motifs(expected)
        actual = motifs.randomized_motif_search(dna, k)

        # compare only score as they exist two different set of motifs with same score
        self.assertEqual(expected_score, actual[1])

    def test_gibs_sampling(self):
        dna = ['CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA',
               'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
               'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
               'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
               'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
        k = 8
        expected = ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TTCAGGTG', 'TCCACGTG']
        expected_score = motifs.score_motifs(expected)
        actual = motifs.gibbs_sampler(dna, k, N=100, times=20)
        self.assertEqual(expected_score, actual[1])

    def test_motif_enumeration_no_mismatch(self):
        dnas = ['CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA',
                'GGGCGAGGTTCTCGGGGGTGCCAAGGTGCCAG',
                'TCTCGGGGAGACCGAAAGAAGTATACAGGCGT',
                'TAGATCAAGTTTCTCGGGGACGTCGGTGAACC',
                'AATCCACCAGCTCCACGTGCAATGTCTCGGGG']
        k = 8
        expected_motifs = {'TCTCGGGG'}
        actual_motifs = motifs.motif_enumeration(dnas, k, d=0)
        self.assertEqual(expected_motifs, actual_motifs)

    def test_motif_enumeration_one_mismatch(self):
        dnas = ['CGCCCCTCTCGAGGGTGTTCAGTAACCGGCCA',
                'GGGCGAGGTTCTCTGGGGTGCCAAGGTGCCAG',
                'TCTCCGGGAGACCGAAAGAAGTATACAGGCGT',
                'TAGATCAAGTTTCTCGGGGACGTCGGTGAACC',
                'AATCCACCAGCTCCACGTGCAATGTCTCAGGG']
        k = 8
        expected_motifs = {'TCTCGGGG'}
        actual_motifs = motifs.motif_enumeration(dnas, k, d=1)
        self.assertEqual(expected_motifs, actual_motifs)

    if __name__ == '__main__':
        unittest.main()
