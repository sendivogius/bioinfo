import functools
import itertools

import math

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


def read_dict_file(name):
    str_dict = dict([line.split() for line in open(name)])
    return {k: int(v) for (k, v) in str_dict.items()}


aminoacids_mass_map = read_dict_file('data/integer_mass_table.txt')
aminoacids_masses = list(sorted(set(aminoacids_mass_map.values())))  # todo why must be sorted?
aminoacids = set(rna_codons.values()) - {''}


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
    aa_len = 3 * len(aa)
    encodings = [rna_to_dna(kmer) for kmer in motifs.get_all_kmers(rna, aa_len, True) if translate(kmer) == aa]
    encodings_rev = [motifs.reverse_strand(rna_to_dna(kmer)) for kmer in
                     motifs.get_all_kmers(rna_rev_comp, aa_len, True) if translate(kmer) == aa]
    return encodings + encodings_rev


def get_peptide_mass(aa_str):
    return sum(aminoacids_mass_map[a] for a in aa_str)


def theoretical_spectrum(aa_seq):
    l = len(aa_seq)
    da = aa_seq * 2
    peptides = [da[start:start + k] for k in range(1, l) for start in range(l)] + ['', aa_seq]
    peptides_mass = [get_peptide_mass(pep) for pep in peptides]
    return sorted(peptides_mass)


def parent_mass(spectrum):
    return max(spectrum)


def expand_peptides(peptides):
    return {p + a for p in peptides for a in aminoacids}


def get_all_combinations_iterative(elements, target_value):
    '''
    Generates ... which sum to target
    :param elements:
    :param target_value:
    :return:
    '''
    element_to_check = elements[0]
    remaining_elements = elements[1:]
    max_num_of_elem = target_value // element_to_check

    if max_num_of_elem * element_to_check == target_value:
        yield {element_to_check: max_num_of_elem}

    if remaining_elements:
        for i in range(max_num_of_elem):
            for current_sol in get_all_combinations_iterative(remaining_elements, target_value - i * element_to_check):
                if i:
                    current_sol[element_to_check] = i
                yield current_sol


def _num_peptide_combinations(peptide_composition):
    pep_len = sum(peptide_composition.values())
    cc = functools.reduce(lambda x, y: x * y, (math.factorial(d) for d in peptide_composition.values()))
    return int(math.factorial(pep_len) / cc)


def number_of_peptides_with_mass(peptide_length):
    comb = get_all_combinations_iterative(aminoacids_masses, peptide_length)
    return sum(_num_peptide_combinations(x) for x in comb)


def is_consistent(test_spec, target_spec):
    return set(test_spec) <= set(target_spec)


def cyclopeptide_sequencing_generator(spectrum):
    peptides = {''}
    while peptides:
        matching_pep = set()
        for pep in expand_peptides(peptides):
            if get_peptide_mass(pep) == parent_mass(spectrum) and theoretical_spectrum(pep) == spectrum:
                yield pep
            if is_consistent(theoretical_spectrum(pep), spectrum):
                matching_pep.add(pep)
        peptides = matching_pep


def mass_string(peptide):
    return '-'.join([str(aminoacids_mass_map[aa]) for aa in peptide])


def cyclopeptide_sequencing(spectrum, unique=False):
    peptides = set(cyclopeptide_sequencing_generator(spectrum))
    if not unique:
        return peptides
    return {mass_string(pep) for pep in peptides}



def approx_pept_mass(start=100, stop=1000, step=2):
    from scipy.optimize import curve_fit
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    import matplotlib.pyplot as plt

    def approx_func(m, k, C):
        return k * np.power(C, m)

    series = pd.Series()
    for i in tqdm(range(start, stop + 1, step)):
        num = number_of_peptides_with_mass(i)
        series.set_value(i, num)
    series.to_csv('peptide_number_by_mass_approx.csv')

    xdata = series.index.values
    ydata = series.values
    popt, pcov = curve_fit(approx_func, xdata, ydata, (0.00742696, 1.02820358))
    print(popt)
    plt.semilogy(xdata, ydata, 'b-', label='data')
    plt.semilogy(xdata, approx_func(xdata, *popt), 'r-', label='approximation')
    plt.xlabel('mass')
    plt.ylabel('number of peptides')
    plt.savefig('a.png')
    plt.show()


if __name__ == "__main__":
    approx_pept_mass(200, 1500, 10)
