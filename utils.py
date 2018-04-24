import collections
from collections import defaultdict
from datetime import datetime
from itertools import takewhile
import operator
from math import log2

import itertools


def reverse_strand(text):
    mappings = {'C': 'G', 'G': 'C', 'A': 'T', 'T': 'A'}
    return ''.join([mappings[nn] for nn in text[::-1]])


def get_all_kmers(pattern, k):
    return {pattern[i:i + k] for i in range(len(pattern) - k + 1)}


def find_occurrences(text, pattern, k=0):
    """get indexes of (possibly overlapping) occurrences of pattern in text"""
    return [i for i in range(len(text) - len(pattern) + 1) if hamming(text[i:i + len(pattern)], pattern) <= k]


def count_occurrences(text, pattern, k=0):
    """count number of (possibly overlapping) occurrences of pattern in text"""
    return len(find_occurrences(text, pattern, k))


def frequent_kmers(text, k, d=0, reverse=False):
    """
    Get most frequent kmers in string with up to d mismatches
    Returns ( {kmer1, kmer2}, num_of_occurences )
    """
    max_kmers = []
    max_kmers_cnt = 0

    reversed = reverse_strand(text) if reverse else ''

    kmers = get_all_kmers(text, k)
    print('kmers', len(kmers))
    kmers_neghbourhood = {p for kmer in kmers for p in get_neighbourhs(kmer, d)}
    print('kmers_neghbourhood', len(kmers_neghbourhood))
    for i, kmer in enumerate(kmers_neghbourhood):
        kmer_counts = count_occurrences(text, kmer, d) + count_occurrences(reversed, kmer, d)
        if kmer_counts == max_kmers_cnt:
            max_kmers.append(kmer)
        elif kmer_counts > max_kmers_cnt:
            max_kmers = [kmer]
            max_kmers_cnt = kmer_counts
        if i % 5000 == 0:
            print(i)
    return set(max_kmers), max_kmers_cnt


def _get_neighbourhs(pattern):
    bases = 'ACTG'
    result = set()
    for i in range(len(pattern)):
        for base in bases:
            result.add(pattern[:i] + base + pattern[(i + 1):])
    return result


def get_neighbourhs(pattern, d):
    result = set([pattern])
    for i in range(d):
        addded = set()
        for pattern in result:
            addded |= _get_neighbourhs(pattern)
        result |= addded
    return result


def get_keys_(dict_, t):
    return {k for (k, v) in dict_.items() if v >= t}


def find_clumps(text, k, L, t):
    """
    Find all distinct k-mers forming (L, t)-clumps in Genome.

    :param text: DNA
    :param k: length of kmers
    :param L: ori length
    :param t: minimum number of occurences of kmer in L
    :return:
    """
    counts = defaultdict(int)
    for k_start in range(L - k + 1):
        counts[text[k_start:k_start + k]] += 1
    kmers = get_keys_(counts, t)

    for L_start in range(1, len(text) - L + 1):
        counts[text[L_start - 1:L_start + k - 1]] -= 1
        new_kmer = text[L_start + L - k:L_start + L]
        counts[new_kmer] += 1
        if counts[new_kmer] >= t:
            kmers.add(new_kmer)

    return kmers


def pattern_to_number(pattern):
    idx = 'ACGT'
    return sum(idx.index(dna_base) * 4 ** i for i, dna_base in enumerate(pattern[::-1]))


def number_to_pattern(number, length):
    idx = 'ACGT'
    pattern = ''
    while number > 0:
        pattern += idx[number % 4]
        number //= 4
    return idx[0] * (length - len(pattern)) + pattern[::-1]


def compute_frequencies(text, k):
    freq_array = [0] * 4 ** k
    for start in range(len(text) - k + 1):
        num = pattern_to_number(text[start:start + k])
        freq_array[num] += 1
    return freq_array


def get_skew(dna, down='C', up='G'):
    skew = [0] * (len(dna) + 1)
    for i, b in enumerate(dna, 1):
        if b == down:
            skew[i] = skew[i - 1] - 1
        elif b == up:
            skew[i] = skew[i - 1] + 1
        else:
            skew[i] = skew[i - 1]
    return skew


def get_min_skew_posiion(dna, down='C', up='G'):
    skew = get_skew(dna, down, up)
    min_skew = min(skew)
    return [pos for (pos, sk) in enumerate(skew) if sk == min_skew]


def plot_skew(dna, down='C', up='G'):
    skew = get_skew(dna, down, up)
    import matplotlib.pyplot as plt
    plt.plot(skew, '-')
    plt.title('Skew plot {}-{}'.format(up, down))
    plt.ylabel('skew')
    plt.xlabel('position')
    plt.show()


def hamming(dna1, dna2):
    assert (len(dna1) == len(dna2))
    return sum(a != b for a, b in zip(dna1, dna2))


def get_profile(motifs, relative=True):
    zeroes = {'A': 0, 'G': 0, 'T': 0, 'C': 0}
    absolute = [{**zeroes, **dict(collections.Counter(m))} for m in zip(*motifs)]
    if relative:
        return [{k: v / total for total in (sum(a.values(), 0.0),) for k, v in a.items()} for a in absolute]
    return absolute


def get_consensus_string(motifs, profile=None):
    profile = profile or get_profile(motifs)
    return ''.join(max(m.items(), key=operator.itemgetter(1))[0] for m in profile)


def score_motifs(motifs, entropy=False):
    if entropy:
        def _calc_entropy(values):
            return sum(-v * log2(v) if v else 0.0 for v in values)

        profile = get_profile(motifs, relative=True)
        return sum(_calc_entropy(val.values()) for val in profile)

    profile = get_profile(motifs, relative=False)
    consensus = get_consensus_string(motifs)
    return sum((sum({**p, **{c: 0}}.values()) for c, p in zip(consensus, profile)))


def dist(pattern, dna):
    def dist_to_single_dna(p, d):
        return min(hamming(p, kmer) for kmer in get_all_kmers(d, len(p)))

    if type(dna) == str:
        dna = [dna]
    return sum(dist_to_single_dna(pattern, d) for d in dna)


def median_string(dna, k):
    all_kmers = map(''.join, itertools.product('ACTG', repeat=k))
    return min(((kmer, dist(kmer, dna)) for kmer in all_kmers), key=operator.itemgetter(1))


if __name__ == "__main__":
    ecoli = open('data\\E_coli.txt').readline()
    plot_skew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')
