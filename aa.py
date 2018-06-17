import itertools

import motifs

rna_codons = {'AAA': 'K',
              'AAC': 'N',
              'AAG': 'K',
              'AAU': 'N',
              'ACA': 'T',
              'ACC': 'T',
              'ACG': 'T',
              'ACU': 'T',
              'AGA': 'R',
              'AGC': 'S',
              'AGG': 'R',
              'AGU': 'S',
              'AUA': 'I',
              'AUC': 'I',
              'AUG': 'M',
              'AUU': 'I',
              'CAA': 'Q',
              'CAC': 'H',
              'CAG': 'Q',
              'CAU': 'H',
              'CCA': 'P',
              'CCC': 'P',
              'CCG': 'P',
              'CCU': 'P',
              'CGA': 'R',
              'CGC': 'R',
              'CGG': 'R',
              'CGU': 'R',
              'CUA': 'L',
              'CUC': 'L',
              'CUG': 'L',
              'CUU': 'L',
              'GAA': 'E',
              'GAC': 'D',
              'GAG': 'E',
              'GAU': 'D',
              'GCA': 'A',
              'GCC': 'A',
              'GCG': 'A',
              'GCU': 'A',
              'GGA': 'G',
              'GGC': 'G',
              'GGG': 'G',
              'GGU': 'G',
              'GUA': 'V',
              'GUC': 'V',
              'GUG': 'V',
              'GUU': 'V',
              'UAA': '',
              'UAC': 'Y',
              'UAG': '',
              'UAU': 'Y',
              'UCA': 'S',
              'UCC': 'S',
              'UCG': 'S',
              'UCU': 'S',
              'UGA': '',
              'UGC': 'C',
              'UGG': 'W',
              'UGU': 'C',
              'UUA': 'L',
              'UUC': 'F',
              'UUG': 'L',
              'UUU': 'F'}

def dna_to_rna(dna):
    return dna.replace('T', 'U')

def rna_to_dna(rna):
    return rna.replace('U', 'T')

def _reverse_dict(d):
    rev = dict()
    for k, v in d.items():
        if v in rev:
            rev[v] += [k]
        else:
            rev[v] = [k]
    return rev


rna_reverse_codons = _reverse_dict(rna_codons)


def translate(dna):
    assert (len(dna) % 3 == 0)
    rna = dna.replace('T', 'U')
    cod = [rna[i:i + 3] for i in range(0, len(rna), 3)]
    aa = [rna_codons[c] for c in cod]
    return ''.join(aa)


def reverse_translate(aa):
    codons_rna = (rna_reverse_codons[a] for a in aa)
    combined = itertools.product(*codons_rna)
    rnas = {''.join(l) for l in combined}
    return rnas


def get_encodings(dna, aa):
    rna = dna_to_rna(dna)
    rna_rev_comp = dna_to_rna(motifs.reverse_strand(dna))
    aa_len = 3*len(aa)
    encodings = [rna_to_dna(kmer) for kmer in motifs.get_all_kmers(rna, aa_len, True) if translate(kmer) == aa]
    encodings_rev = [motifs.reverse_strand(rna_to_dna(kmer)) for kmer in motifs.get_all_kmers(rna_rev_comp, aa_len, True) if translate(kmer) == aa]
    return encodings + encodings_rev


if __name__ == "__main__":
    bacillus = ''.join((s.strip() for s in open('data\\Bacillus_brevis.txt').readlines()))
    tyrocidineB1 = 'VKLFPWFNQY'
    encodings = get_encodings(bacillus, tyrocidineB1)
    print(encodings)
