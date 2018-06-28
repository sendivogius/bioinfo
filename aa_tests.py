import unittest
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

    def test_peptide_mass(self):
        aa_str = 'NQE'
        expected = 371
        self.assertEqual(expected, aa.get_peptide_mass(aa_str))

    def test_peptide_mass_empty(self):
        aa_str = ''
        expected = 0
        self.assertEqual(expected, aa.get_peptide_mass(aa_str))

    def test_aa_spectrum(self):
        aa_str = 'LEQN'
        expected = [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
        self.assertEqual(expected, aa.theoretical_spectrum(aa_str))

    def test_parent_mass(self):
        spectrum = [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
        expected = 484
        self.assertEqual(expected, aa.parent_mass(spectrum))

    def test_expand_peptides_empty(self):
        peptides = {''}
        expected = {"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"}
        self.assertEqual(expected, aa.expand_peptides(peptides))

    def test_expand_peptides(self):
        peptides = {'A', 'DD'}
        expected = {"AA", "AR", "AN", "AD", "AC", "AE", "AQ", "AG", "AH", "AI", "AL", "AK", "AM", "AF", "AP", "AS",
                    "AT", "AW", "AY", "AV",
                    "DDA", "DDR", "DDN", "DDD", "DDC", "DDE", "DDQ", "DDG", "DDH", "DDI", "DDL", "DDK", "DDM", "DDF",
                    "DDP", "DDS", "DDT", "DDW", "DDY", "DDV"}
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

    def test_cyclopeptide_sequencing_unique(self):
        spectrum = [0,71,97,99,103,113,113,114,115,131,137,196,200,202,208,214,226,227,228,240,245,299,311,311,316,327,337,339,340,341,358,408,414,424,429,436,440,442,453,455,471,507,527,537,539,542,551,554,556,566,586,622,638,640,651,653,657,664,669,679,685,735,752,753,754,756,766,777,782,782,794,848,853,865,866,867,879,885,891,893,897,956,962,978,979,980,980,990,994,996,1022,1093]
        print( aa.cyclopeptide_sequencing(spectrum, unique=True))
        self.assertEqual({'113-128-186', '113-186-128', '128-113-186', '128-186-113', '186-113-128', '186-128-113'},
                         aa.cyclopeptide_sequencing(spectrum, unique=True))

    def test_cyclopeptide_sequencing(self):
        spectrum = [11,71,87,97,99,101,115,115,128,128,129,172,184,186,186,198,227,243,243,244,257,269,283,285,287,314,314,342,356,372,372,372,384,384,411,415,429,443,455,471,471,487,500,512,512,526,544,558,558,570,583,599,599,615,627,641,655,659,686,686,698,698,698,714,728,756,756,783,785,787,801,813,826,827,827,843,872,884,884,886,898,941,942,942,955,955,969,971,973,983,999,1070]

        self.assertEqual(
            {'IWQ', 'WQL', 'WKL', 'QWL', 'IQW', 'KLW', 'LWK', 'WQI', 'LWQ', 'WLQ', 'KIW', 'WIQ', 'LKW', 'IKW', 'KWI',
             'WIK', 'LQW', 'QWI', 'IWK', 'KWL', 'QLW', 'WLK', 'WKI', 'QIW'}
            , aa.cyclopeptide_sequencing(spectrum, unique=False))

    def test_a(self):
        pep = 'CHAMNLLDVP'
        spec = aa.theoretical_spectrum(pep)

        target = [0,71,97,99,103,113,113,114,115,131,137,196,200,202,208,214,226,227,228,240,245,299,311,311,316,327,337,339,340,341,358,408,414,424,429,436,440,442,453,455,471,507,527,537,539,542,551,554,556,566,586,622,638,640,651,653,657,664,669,679,685,735,752,753,754,756,766,777,782,782,794,848,853,865,866,867,879,885,891,893,897,956,962,978,979,980,980,990,994,996,1022,1093]
        cons = aa.is_consistent(spec, target)
        self.assertTrue(cons)

if __name__ == '__main__':
    unittest.main()
