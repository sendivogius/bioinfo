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
        print(aa.get_encodings(dna, aa_str))
        self.assertEqual(sorted(expected), sorted(aa.get_encodings(dna, aa_str)))

if __name__ == '__main__':
    unittest.main()
