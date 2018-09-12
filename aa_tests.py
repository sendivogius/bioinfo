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
        self.assertEqual(expected, aa.get_parent_mass(spectrum))

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
        spectrum = [0, 113, 128, 186, 241, 299, 314, 427]
        # print(aa.cyclopeptide_sequencing(spectrum, unique=True))
        self.assertEqual({'113-128-186', '113-186-128', '128-113-186', '128-186-113', '186-113-128', '186-128-113'},
                         aa.cyclopeptide_sequencing(spectrum, unique=True))

    def test_cyclopeptide_sequencing(self):
        spectrum = [0, 113, 128, 186, 241, 299, 314, 427]

        self.assertEqual(
            {'IWQ', 'WQL', 'WKL', 'QWL', 'IQW', 'KLW', 'LWK', 'WQI', 'LWQ', 'WLQ', 'KIW', 'WIQ', 'LKW', 'IKW', 'KWI',
             'WIK', 'LQW', 'QWI', 'IWK', 'KWL', 'QLW', 'WLK', 'WKI', 'QIW'},
            aa.cyclopeptide_sequencing(spectrum, unique=False))

    def test_peptide_spectrum_space(self):
        peptide = 'NQEL'
        spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
        self.assertEqual(11, aa.score_peptide(spectrum, peptide))

    def test_peptide_spectrum_space_multipled1(self):
        peptide = 'NNEQ'
        spectrum = [0, 114, 10, 114, 432, 114]
        self.assertEqual(3, aa.score_peptide(spectrum, peptide))

    def test_leadeboard_cyclopeptide_sequencing(self):
        N = 10
        spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]
        N = 300
        spectrum = [0,87,87,87,97,97,99,99,99,99,101,101,103,113,113,114,115,115,115,129,129,131,137,147,147,174,184,196,200,200,204,212,212,214,215,216,218,218,228,228,230,230,230,236,244,246,250,261,276,287,303,305,305,311,313,315,317,318,325,327,329,329,329,331,343,345,349,349,360,362,373,375,383,390,392,402,404,418,426,426,428,431,440,442,442,444,448,450,458,460,461,462,465,472,489,490,491,496,497,505,518,521,525,541,541,541,541,549,549,555,559,564,573,578,579,587,591,595,598,602,605,605,608,610,622,636,638,640,640,656,665,670,670,674,677,678,678,678,690,696,701,702,705,709,709,711,720,736,749,752,753,755,764,765,769,769,771,775,777,777,789,805,806,810,810,814,814,817,817,823,825,851,852,864,868,870,870,874,878,883,890,892,898,901,904,911,913,914,916,920,920,924,927,939,953,969,970,977,979,981,982,983,988,989,999,1001,1005,1007,1013,1014,1014,1021,1026,1035,1038,1040,1045,1067,1067,1069,1080,1082,1082,1092,1094,1098,1100,1101,1104,1113,1119,1120,1122,1127,1132,1135,1136,1139,1142,1143,1146,1166,1181,1195,1196,1200,1206,1207,1209,1211,1214,1214,1219,1219,1219,1229,1231,1232,1236,1242,1243,1245,1250,1251,1256,1295,1296,1301,1303,1306,1310,1313,1318,1318,1319,1322,1330,1331,1339,1342,1343,1343,1348,1351,1355,1358,1360,1365,1366,1405,1410,1411,1416,1418,1419,1425,1429,1430,1432,1442,1442,1442,1447,1447,1450,1452,1454,1455,1461,1465,1466,1480,1495,1515,1518,1519,1522,1525,1526,1529,1534,1539,1541,1542,1548,1557,1560,1561,1563,1567,1569,1579,1579,1581,1592,1594,1594,1616,1621,1623,1626,1635,1640,1647,1647,1648,1654,1656,1660,1662,1672,1673,1678,1679,1680,1682,1684,1691,1692,1708,1722,1734,1737,1741,1741,1745,1747,1748,1750,1757,1760,1763,1769,1771,1778,1783,1787,1791,1791,1793,1797,1809,1810,1836,1838,1844,1844,1847,1847,1851,1851,1855,1856,1872,1884,1884,1886,1890,1892,1892,1896,1897,1906,1908,1909,1912,1925,1941,1950,1952,1952,1956,1959,1960,1965,1971,1983,1983,1983,1984,1987,1991,1991,1996,2005,2021,2021,2023,2025,2039,2051,2053,2056,2056,2059,2063,2066,2070,2074,2082,2083,2088,2097,2102,2106,2112,2112,2120,2120,2120,2120,2136,2140,2143,2156,2164,2165,2170,2171,2172,2189,2196,2199,2200,2201,2203,2211,2213,2217,2219,2219,2221,2230,2233,2235,2235,2243,2257,2259,2269,2271,2278,2286,2288,2299,2301,2312,2312,2318,2330,2332,2332,2332,2334,2336,2343,2344,2346,2348,2350,2356,2356,2358,2374,2385,2400,2411,2415,2417,2425,2431,2431,2431,2433,2433,2443,2443,2445,2446,2447,2449,2449,2457,2461,2461,2465,2477,2487,2514,2514,2524,2530,2532,2532,2546,2546,2546,2547,2548,2548,2558,2560,2560,2562,2562,2562,2562,2564,2564,2574,2574,2574,2661]
        self.assertEqual('FAEL', aa.leaderboard_cyclopeptide_sequencing(spectrum, N, True))

    def test_trim_leaderboard(self):
        spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]
        leaderboard = {'A', 'S', 'AF'}
        N = 2
        self.assertEqual({'A', 'AF'}, aa.trim_leaderboard(leaderboard, N, spectrum))

    def test_trim_leaderboard_tie(self):
        spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]
        leaderboard = {'A', 'F', 'S', 'AF'}
        N = 2
        self.assertEqual({'A', 'AF', 'F'}, aa.trim_leaderboard(leaderboard, N, spectrum))

    def test_get_n_for_ties_1(self):
        elems = [1, 2, 3]
        n = 2
        self.assertEqual(2, aa._get_n_for_ties(elems, n,))

    def test_get_n_for_ties_2(self):
        elems = [1, 2, 2, 3]
        n = 2
        self.assertEqual(3, aa._get_n_for_ties(elems, n, ))

    def test_get_n_for_ties_3(self):
        elems = [1, 2, 2, 2]
        n = 2
        self.assertEqual(4, aa._get_n_for_ties(elems, n, ))

if __name__ == '__main__':
    unittest.main()
