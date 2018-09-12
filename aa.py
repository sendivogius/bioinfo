import functools
import itertools

import math
from collections import Counter
from operator import itemgetter

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
aminoacids_mass_reverse_map = _reverse_dict(aminoacids_mass_map)


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
    if type(aa_str) == str:
        return sum(aminoacids_mass_map[a] for a in aa_str)
    return sum(aa_str)


def theoretical_spectrum(aa_seq, cyclic=True):
    "Generates cyclospectrum of peptide"
    l = len(aa_seq)
    if type(aa_seq) == str:
        empty = ''
    else:
        empty = ()
    da = aa_seq
    if cyclic:
        da = aa_seq * 2
    peptides = [da[start:start + k] for k in range(1, l) for start in range(l) if cyclic or (start + k) <= l] + [empty,
                                                                                                                 aa_seq]
    peptides_mass = (get_peptide_mass(pep) for pep in peptides)
    return Counter(peptides_mass)


def get_parent_mass(spectrum):
    return max(spectrum)


def expand_peptides(mass_tuples):
    return {p + (m,) for p in mass_tuples for m in aminoacids_masses}


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
    parent_mass = get_parent_mass(spectrum)
    peptides = {()}
    while peptides:
        matching_pep = set()
        for pep in expand_peptides(peptides):
            if get_peptide_mass(pep) == parent_mass and theoretical_spectrum(pep) == spectrum:
                yield pep
            if is_consistent(theoretical_spectrum(pep), spectrum):
                matching_pep.add(pep)
        peptides = matching_pep


def mass_string(peptide):
    return '-'.join([str(aminoacids_mass_map[aa]) for aa in peptide])


def cyclopeptide_sequencing(spectrum, integers=False):
    peptides = set(cyclopeptide_sequencing_generator(spectrum))
    if integers:
        return peptides
    a = itertools.chain.from_iterable(masses_to_peptide(pep) for pep in peptides)
    return set(a)


def _get_n_for_ties(elements, N):
    len_elem = len(elements)
    if len_elem <= N:
        return len_elem
    new_n = N
    while new_n < len_elem and elements[new_n] == elements[N - 1]:
        new_n += 1
    return new_n


def trim_leaderboard(leaderboard, N, spectrum):
    if not leaderboard:
        return set(), 0, 0
    cyclic_score = N == 1  # if N == 1 we ge best which should be scored diffrently
    scores_tuples = sorted([(score_peptide(spectrum, l, cyclic=False), l) for l in leaderboard], reverse=True)
    sorted_scores = list(map(itemgetter(0), scores_tuples))
    new_n = _get_n_for_ties(sorted_scores, N)
    return set(s[1] for s in scores_tuples[:new_n]), sorted_scores[0], sorted_scores[new_n - 1]


def leaderboard_cyclopeptide_sequencing(spectrum, N=50, unique=False):
    parent_mass = get_parent_mass(spectrum)
    leaderboard = {''}
    leader_peptides = ''
    leader_score = 0
    while leaderboard:
        new_board = set()
        for pep in expand_peptides(leaderboard):
            pep_mass = get_peptide_mass(pep)
            if pep_mass == parent_mass:
                score = score_peptide(spectrum, pep, cyclic=True)
                if score > leader_score:
                    leader_peptides = {pep}
                    leader_score = score
                elif score == leader_score:
                    leader_peptides.add(pep)
            if pep_mass <= parent_mass:
                new_board.add(pep)
        leaderboard, _, _ = trim_leaderboard(new_board, N, spectrum)

    if not unique:
        return leader_peptides, leader_score
    return {mass_string(pep) for pep in leader_peptides}, leader_score


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

def masses_to_peptide(mass_tuple):
    aminoacids = (aminoacids_mass_reverse_map[m] for m in mass_tuple)
    combined = itertools.product(*aminoacids)
    peptides = {''.join(l) for l in combined}
    return peptides

def score_peptide(spectrum, peptide, cyclic=True):
    # print(peptide)
    peptide_spectrum = theoretical_spectrum(peptide, cyclic)
    common_masses = spectrum.keys() & peptide_spectrum.keys()
    score = sum(min(spectrum[k], peptide_spectrum[k]) for k in common_masses)
    return score


if __name__ == "__main__":
    # approx_pept_mass(200, 1500, 10)
    pep = [97, 71, 115, 147, 114, 128, 163, 99, 128, 113, 147]
    import pprint

    spec_10 = Counter(
        [0, 97, 99, 114, 128, 147, 147, 163, 186, 227, 241, 242, 244, 260, 261, 262, 283, 291, 333, 340, 357, 385,
         389, 390, 390, 405, 430, 430, 447, 485, 487, 503, 504, 518, 543, 544, 552, 575, 577, 584, 632, 650, 651,
         671, 672, 690, 691, 738, 745, 747, 770, 778, 779, 804, 818, 819, 820, 835, 837, 875, 892, 917, 932, 932,
         933, 934, 965, 982, 989, 1030, 1039, 1060, 1061, 1062, 1078, 1080, 1081, 1095, 1136, 1159, 1175, 1175, 1194,
         1194, 1208, 1209, 1223, 1225, 1322])
    spec_25 = Counter(
        [1, 309, 330, 333, 340, 347, 385, 388, 389, 390, 390, 405, 435, 447, 485, 487, 503, 504, 518, 544, 552, 575,
         577, 584, 599, 608, 631, 632, 650, 651, 653, 672, 690, 691, 717, 738, 745, 770, 779, 804, 818, 819, 827, 835,
         837, 875, 892, 892, 917, 932, 932, 933, 934, 965, 982, 989, 1039, 1060, 1062, 1078, 1080, 1081, 1095, 1136,
         1159, 1175, 1175, 1194, 1194, 1208, 1209, 1223, 1322])
    N = 1000

    # best_10 = leaderboard_cyclopeptide_sequencing(spec_10, N, True)
    # pprint.pprint(best_10)
    best_25 = leaderboard_cyclopeptide_sequencing(spec_25, N, True)
    pprint.pprint(best_25)
