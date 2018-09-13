import unittest
from collections import Counter

import aa


class AaTests(unittest.TestCase):
    def test_translation(self):
        dna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
        expected = 'MAMAPRTEINSTRING'
        self.assertEqual(expected, aa.translate(dna))

    def test_reverse_translation(self):
        aa_str = 'MS'
        expected = {
            'AUGAGU',
            'AUGAGC',
            'AUGUCA',
            'AUGUCC',
            'AUGUCG',
            'AUGUCU'
        }
        self.assertEqual(expected, aa.reverse_translate(aa_str))

    def test_get_econdings(self):
        dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
        aa_str = 'MA'
        expected = [
            'ATGGCC',
            'GGCCAT',
            'ATGGCC',
        ]
        self.assertEqual(sorted(expected), sorted(aa.get_encodings(dna, aa_str)))

    def test_peptide_mass_str(self):
        aa_str = 'NQE'
        expected = 371
        self.assertEqual(expected, aa.get_peptide_mass(aa_str))

    def test_peptide_mass_empty_str(self):
        aa_str = ''
        expected = 0
        self.assertEqual(expected, aa.get_peptide_mass(aa_str))

    def test_peptide_mass_tuple(self):
        mass_tuple = (128, 129, 114)
        expected = 371
        self.assertEqual(expected, aa.get_peptide_mass(mass_tuple))

    def test_peptide_mass_empty(self):
        mass_tuple = ()
        expected = 0
        self.assertEqual(expected, aa.get_peptide_mass(mass_tuple))

    def test_aa_spectrum_cyclic_str(self):
        aa_str = 'LEQN'
        expected = Counter([0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484])
        self.assertEqual(expected, aa.theoretical_spectrum(aa_str, cyclic=True))

    def test_aa_spectrum_linear_str(self):
        aa_str = 'LEQN'
        expected = Counter([113, 129, 128, 114, 242, 257, 242, 370, 371, 0, 484])
        self.assertEqual(expected, aa.theoretical_spectrum(aa_str, cyclic=False))

    def test_aa_spectrum_cyclic_tuple(self):
        mass_tuple = (113, 129, 128, 114)
        expected = Counter([0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484])
        self.assertEqual(expected, aa.theoretical_spectrum(mass_tuple, cyclic=True))

    def test_aa_spectrum_linear_tuple(self):
        mass_tuple = (113, 129, 128, 114)
        expected = Counter([113, 129, 128, 114, 242, 257, 242, 370, 371, 0, 484])
        self.assertEqual(expected, aa.theoretical_spectrum(mass_tuple, cyclic=False))

    def test_parent_mass(self):
        spectrum = Counter([0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484])
        expected = 484
        self.assertEqual(expected, aa.get_parent_mass(spectrum))

    def test_expand_peptides_empty(self):
        peptides = {()}
        expected = {(57,), (71,), (87,), (97,), (99,), (101,), (103,), (113,), (114,), (115,), (128,), (129,), (131,),
                    (137,), (147,), (156,), (163,), (186,), }
        self.assertEqual(expected, aa.expand_peptides(peptides))

    def test_expand_peptides(self):
        peptides = {(71,), (115, 115)}
        expected = {(71, 57), (71, 71), (71, 87), (71, 97), (71, 99), (71, 101), (71, 103), (71, 113), (71, 114),
                    (71, 115), (71, 128), (71, 129), (71, 131), (71, 137), (71, 147), (71, 156), (71, 163), (71, 186),
                    (115, 115, 57), (115, 115, 71), (115, 115, 87), (115, 115, 97), (115, 115, 99), (115, 115, 101),
                    (115, 115, 103), (115, 115, 113), (115, 115, 114), (115, 115, 115), (115, 115, 128),
                    (115, 115, 129), (115, 115, 131), (115, 115, 137), (115, 115, 147), (115, 115, 156),
                    (115, 115, 186), (115, 115, 163), }
        self.assertEqual(expected, aa.expand_peptides(peptides))

    def test_subpeptides_number(self):
        peptide_length = 841
        expected = 112663579
        self.assertEqual(expected, aa.number_of_peptides_with_mass(peptide_length))

    def test_is_not_consistent(self):
        test_peptide = 'VKF'
        test_spec = aa.theoretical_spectrum(test_peptide)
        target_spec = [0, 97, 99, 113, 114, 128, 128, 147, 147, 163, 186, 227, 241, 242, 244, 260, 261, 262, 283, 291,
                       333, 340, 357, 388, 389, 390, 390, 405, 430, 430, 447, 485, 487, 503, 504, 518, 543, 544, 552,
                       575, 577, 584, 631, 632, 650, 651, 671, 672, 690, 691, 738, 745, 747, 770, 778, 779, 804, 818,
                       819, 835, 837, 875, 892, 892, 917, 932, 932, 933, 934, 965, 982, 989, 1031, 1039, 1060, 1061,
                       1062, 1078, 1080, 1081, 1095, 1136, 1159, 1175, 1175, 1194, 1194, 1208, 1209, 1223, 1225, 1322]
        self.assertFalse(aa.is_consistent(test_spec, target_spec))

    def test_is_consistent(self):
        test_peptide = 'VKY'
        test_spec = aa.theoretical_spectrum(test_peptide)
        target_spec = [0, 97, 99, 113, 114, 128, 128, 147, 147, 163, 186, 227, 241, 242, 244, 260, 261, 262, 283, 291,
                       333, 340, 357, 388, 389, 390, 390, 405, 430, 430, 447, 485, 487, 503, 504, 518, 543, 544, 552,
                       575, 577, 584, 631, 632, 650, 651, 671, 672, 690, 691, 738, 745, 747, 770, 778, 779, 804, 818,
                       819, 835, 837, 875, 892, 892, 917, 932, 932, 933, 934, 965, 982, 989, 1031, 1039, 1060, 1061,
                       1062, 1078, 1080, 1081, 1095, 1136, 1159, 1175, 1175, 1194, 1194, 1208, 1209, 1223, 1225, 1322]
        self.assertTrue(aa.is_consistent(test_spec, target_spec))

    def test_mass_string(self):
        peptide = 'IWK'
        expected = '128-111-111'
        self.assertTrue(expected, aa.mass_string(peptide))

    def test_cyclopeptide_sequencing_integers(self):
        spectrum = Counter([0, 113, 128, 186, 241, 299, 314, 427])
        self.assertEqual(
            {(113, 128, 186), (113, 186, 128), (128, 113, 186), (128, 186, 113), (186, 113, 128), (186, 128, 113)},
            aa.cyclopeptide_sequencing(spectrum, integers=True))

    def test_cyclopeptide_sequencing(self):
        spectrum = Counter([0, 113, 128, 186, 241, 299, 314, 427])
        self.assertEqual(
            {'IWQ', 'WQL', 'WKL', 'QWL', 'IQW', 'KLW', 'LWK', 'WQI', 'LWQ', 'WLQ', 'KIW', 'WIQ', 'LKW', 'IKW', 'KWI',
             'WIK', 'LQW', 'QWI', 'IWK', 'KWL', 'QLW', 'WLK', 'WKI', 'QIW'},
            aa.cyclopeptide_sequencing(spectrum, integers=False))

    def test_peptide_score_cyclic(self):
        peptide = 'NQEL'
        spectrum = Counter([0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484])
        self.assertEqual(11, aa.score_peptide(spectrum, peptide, cyclic=True))

    def test_peptide_score_linear(self):
        peptide = 'NQEL'
        spectrum = Counter([0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484])
        self.assertEqual(8, aa.score_peptide(spectrum, peptide, cyclic=False))

    def test_peptide_peptide_score_multipled1(self):
        peptide = 'NNEQ'
        spectrum = Counter([0, 114, 10, 114, 432, 114])
        self.assertEqual(3, aa.score_peptide(spectrum, peptide))

    def test_leadeboard_cyclopeptide_sequencing_integers(self):
        N = 10
        spectrum = Counter([0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460])
        self.assertEqual(({(71, 129, 113, 147),
                           (71, 147, 113, 129),
                           (113, 129, 71, 147),
                           (113, 147, 71, 129),
                           (129, 71, 147, 113),
                           (129, 113, 147, 71),
                           (147, 71, 129, 113),
                           (147, 113, 129, 71)},
                          13),
                         aa.leaderboard_cyclopeptide_sequencing(spectrum, N, integers=True))

    def test_leadeboard_cyclopeptide_sequencing(self):
        N = 10
        spectrum = Counter([0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460])
        self.assertEqual(({'AEIF', 'AELF', 'AFIE', 'AFLE', 'EAFI', 'EAFL', 'EIFA', 'ELFA', 'FAEI', 'FAEL', 'FIEA',
                           'FLEA', 'IEAF', 'IFAE', 'LEAF', 'LFAE'},
                          13),
                         aa.leaderboard_cyclopeptide_sequencing(spectrum, N, integers=False))

    def test_trim_leaderboard(self):
        spectrum = Counter([0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460])
        leaderboard = {'A', 'S', 'AF'}
        N = 2
        self.assertEqual(({'A', 'AF'}, 4, 2), aa.trim_leaderboard(leaderboard, N, spectrum))

    def test_trim_leaderboard_tie(self):
        spectrum = Counter([0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460])
        leaderboard = {'A', 'F', 'S', 'AF'}
        N = 2
        self.assertEqual(({'A', 'AF', 'F'}, 4, 2), aa.trim_leaderboard(leaderboard, N, spectrum))

    def test_get_n_for_ties_1(self):
        elems = [1, 2, 3]
        n = 2
        self.assertEqual(2, aa._get_n_for_ties(elems, n, ))

    def test_get_n_for_ties_2(self):
        elems = [1, 2, 2, 3]
        n = 2
        self.assertEqual(3, aa._get_n_for_ties(elems, n, ))

    def test_get_n_for_ties_3(self):
        elems = [1, 2, 2, 2]
        n = 2
        self.assertEqual(4, aa._get_n_for_ties(elems, n, ))

    def test_masses_to_peptide(self):
        spectrum = (97, 113, 147, 128)
        expected = {'PIFK', 'PLFK', 'PIFQ', 'PLFQ'}
        self.assertEqual(expected, aa.masses_to_peptide(spectrum))


if __name__ == '__main__':
    unittest.main()
