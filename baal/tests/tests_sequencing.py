import unittest
from baal import sequencing


class SequencingTests(unittest.TestCase):
    def test_composition(self):
        pattern = 'TATGGGGTGC'
        k = 3
        expected = ['ATG', 'GGG', 'GGG', 'GGT', 'GTG', 'TAT', 'TGC', 'TGG']
        self.assertEqual(expected, sequencing.composition(pattern, k))

    def test_genome_from_sequence(self):
        seq = ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']
        expected = 'ACCGAAGCT'
        self.assertEqual(expected, sequencing.genome_from_sequence(seq))

    def test_genome_from_sequence2(self):
        seq = ['GGC', 'GCT', 'CTT', 'TTA', 'TAC', 'ACC', 'CCA']
        expected = 'GGCTTACCA'
        self.assertEqual(expected, sequencing.genome_from_sequence(seq))

    def test_graph_from_sequence(self):
        seq = ['ATGCG', 'GCATG', 'CATGC', 'AGGCA', 'GGCAT', 'GGCAC']
        expected = {'AGGCA': ['GGCAC', 'GGCAT'],
                    'CATGC': ['ATGCG'],
                    'GCATG': ['CATGC'],
                    'GGCAT': ['GCATG']}
        self.assertEqual(expected, sequencing.graph_from_sequence(seq))

    def test_debrujin_graph_from_sequence(self):
        dna = 'AAGATTCTCTAAGA'
        k = 4
        expected = {'AAG': ['AGA', 'AGA'],
                    'AGA': ['GAT'],
                    'ATT': ['TTC'],
                    'CTA': ['TAA'],
                    'CTC': ['TCT'],
                    'GAT': ['ATT'],
                    'TAA': ['AAG'],
                    'TCT': ['CTA', 'CTC'],
                    'TTC': ['TCT']}
        self.assertEqual(expected, sequencing.debrujin_from_dna(dna, k))

    def test_debrujin_graph_from_kmers(self):
        kmers = ['GAGG',
                 'CAGG',
                 'GGGG',
                 'GGGA',
                 'CAGG',
                 'AGGG',
                 'GGAG']
        expected = {
            'AGG': ['GGG'],
            'CAG': ['AGG', 'AGG'],
            'GAG': ['AGG'],
            'GGA': ['GAG'],
            'GGG': ['GGA', 'GGG']
        }

        self.assertEqual(expected, sequencing.debrujin_from_kmers(kmers))

    def test_eulerian_cycle(self):
        graph = {
            '0': ['3'],
            '1': ['0'],
            '2': ['1', '6'],
            '3': ['2'],
            '4': ['2'],
            '5': ['4'],
            '6': ['5', '8'],
            '7': ['9'],
            '8': ['7'],
            '9': ['6'],

        }
        expected = ['0', '3', '2', '6', '8', '7', '9', '6', '5', '4', '2', '1', '0']
        self.assertEqual(expected, sequencing.eulerian_cycle(graph))

    def test_eulerian_path(self):
        graph = {
            '0': ['2'],
            '1': ['3'],
            '2': ['1'],
            '3': ['4', '0'],
            '6': ['3', '7'],
            '7': ['8'],
            '8': ['9'],
            '9': ['6'],
        }

        expected = ['6', '7', '8', '9', '6', '3', '0', '2', '1', '3', '4']
        self.assertEqual(expected, sequencing.eulerian_path(graph))

    def test_string_reconstruction(self):
        kmers = ['CTTA', 'ACCA', 'TACC', 'GGCT', 'GCTT', 'TTAC']
        expected = 'GGCTTACCA'
        self.assertEqual(expected, sequencing.string_reconstruction(kmers))

    def test_universal_string(self):
        k = 4
        expected = '0001111011001010'
        self.assertEqual(expected, sequencing.universal_string(k))

    def test_paired_composition(self):
        dna = 'TAATGCCATGGGATGTT'
        k = 3
        d = 1
        expected = [('AAT', 'CCA'), ('ATG', 'CAT'), ('ATG', 'GAT'), ('CAT', 'GGA'), ('CCA', 'GGG'), ('GCC', 'TGG'),
                    ('GGA', 'GTT'), ('GGG', 'TGT'), ('TAA', 'GCC'), ('TGC', 'ATG'), ('TGG', 'ATG')]
        self.assertEqual(expected, sequencing.paired_composition(dna, k, d))

    def test_kmer_prefix_str(self):
        kmer = 'ACTG'
        expected = 'ACT'
        self.assertEqual(expected, sequencing.kmer_prefix(kmer))

    def test_kmer_prefix_tuple(self):
        kmers = ('ACTG', 'GGAC')
        expected = ('ACT', 'GGA')
        self.assertEqual(expected, sequencing.kmer_prefix(kmers))

    def test_kmer_suffix_str(self):
        kmer = 'ACTG'
        expected = 'CTG'
        self.assertEqual(expected, sequencing.kmer_suffix(kmer))

    def test_kmer_suffix_tuple(self):
        kmers = ('ACTG', 'GGAC')
        expected = ('CTG', 'GAC')
        self.assertEqual(expected, sequencing.kmer_suffix(kmers))

    def test_debrujin_graph_from_kmers_pairs(self):
        kmers = [
            ('AG', 'AG'),
            ('GC', 'GC'),
            ('CA', 'CT'),
            ('AG', 'TG'),
            ('GC', 'GC'),
            ('CT', 'CT'),
            ('TG', 'TG'),
            ('GC', 'GC'),
            ('CT', 'CA'),
        ]
        expected = {('A', 'A'): [('G', 'G')],
                    ('A', 'T'): [('G', 'G')],
                    ('C', 'C'): [('A', 'T'), ('T', 'A'), ('T', 'T')],
                    ('G', 'G'): [('C', 'C'), ('C', 'C'), ('C', 'C')],
                    ('T', 'T'): [('G', 'G')]
                    }
        self.assertEqual(expected, sequencing.debrujin_from_kmers(kmers))

    def test_string_reconstruction_pair(self):
        kmers = [('GAGA', 'TTGA'),
                 ('TCGT', 'GATG'),
                 ('CGTG', 'ATGT'),
                 ('TGGT', 'TGAG'),
                 ('GTGA', 'TGTT'),
                 ('GTGG', 'GTGA'),
                 ('TGAG', 'GTTG'),
                 ('GGTC', 'GAGA'),
                 ('GTCG', 'AGAT')]
        d = 2
        expected = 'GTGGTCGTGAGATGTTGA'
        self.assertEqual(expected, sequencing.string_pair_reconstruction(kmers, d))

    def test_string_reconstruction_pair_multiple_paths(self):
        kmers = [('AG', 'AG'), ('GC', 'GC'), ('CA', 'CT'), ('AG', 'TG'), ('GC', 'GC'), ('CT', 'CT'), ('TG', 'TG'),
                 ('GC', 'GC'), ('CT', 'CA')]
        d = 1
        expected = 'AGCAGCTGCTGCA'
        self.assertEqual(expected, sequencing.string_pair_reconstruction(kmers, d))

    def test_genome_from_sequence_pair(self):
        seq = [('GTG', 'GTG'),
               ('TGG', 'TGA'),
               ('GGT', 'GAG'),
               ('GTC', 'AGA'),
               ('TCG', 'GAT'),
               ('CGT', 'ATG'),
               ('GTG', 'TGT'),
               ('TGA', 'GTT'),
               ('GAG', 'TTG'),
               ('AGA', 'TGA')]
        d = 3
        expected = 'GTGGTCGTGAGATGTTGA'
        self.assertEqual(expected, sequencing.genome_from_pair_sequence(seq, d))

    def test_genome_from_sequence_pair2(self):
        seq = [('A', 'A'),
               ('G', 'G'),
               ('C', 'C'),
               ('A', 'T'),
               ('G', 'G'),
               ('C', 'C'),
               ('T', 'T'),
               ('G', 'G'),
               ('C', 'C'),
               ('T', 'A')]
        d = 2
        expected = 'AGCAGCTGCTGCA'
        self.assertEqual(expected, sequencing.genome_from_pair_sequence(seq, d))

    def test_genome_from_sequence_pair_invalid(self):
        seq = [('A', 'A'),
               ('G', 'G'),
               ('C', 'C'),
               ('T', 'T'),
               ('G', 'G'),
               ('C', 'C'),
               ('A', 'T'),
               ('G', 'G'),
               ('C', 'C'),
               ('T', 'A')]
        d = 2
        expected = None
        self.assertEqual(expected, sequencing.genome_from_pair_sequence(seq, d))

    def test_get_contigs(self):
        kmers = ['ATG', 'ATG', 'TGT', 'TGG', 'CAT', 'GGA', 'GAT', 'AGA']
        expected = ['ATG', 'ATG', 'TGT', 'TGGA', 'CAT', 'GAT', 'AGA']
        self.assertEqual(sorted(expected), sorted(sequencing.get_contigs(kmers)))

    def test_maximal_nonbranching_paths(self):
        graph = {1: [2],
                 2: [3],
                 3: [4, 5],
                 6: [7],
                 7: [6],
                 10: [11],
                 11: [12],
                 12: [10]
                 }
        expected = [[1, 2, 3], [3, 4], [3, 5], [6, 7, 6], [10, 11, 12, 10]]
        self.assertEqual(sorted(expected), sorted(list((sequencing.maximal_nonbranching_paths(graph)))))


if __name__ == '__main__':
    unittest.main()
