import collections
from collections import defaultdict
from datetime import datetime
from itertools import takewhile
import operator


def reverse_strand(text):
    mappings = {'C': 'G', 'G': 'C', 'A': 'T', 'T': 'A'}
    return ''.join([mappings[nn] for nn in text[::-1]])


def find_occurrences(text, pattern):
    """get indexes of (possibly overlapping) occurrences of pattern in text"""
    return [i for i in range(len(text) - len(pattern) + 1) if text[i:i + len(pattern)] == pattern]


def count_occurrences(text, pattern):
    """count number of (possibly overlapping) occurrences of pattern in text"""
    return len(find_occurrences(text, pattern))


def frequent_kmers(text, k):
    """
    Get most frequent kmers in string
    Returns ( {kmer1, kmer2}, num_of_occurences )
    """
    l = list(text[i:i + k] for i in range(len(text) - k + 1))
    counts = collections.Counter(l)
    max_cnt = counts.most_common(1)[0][1]
    return set(map(operator.itemgetter(0), takewhile(lambda x: x[1] >= max_cnt, counts.most_common()))), max_cnt


def find_clumps(text, k, L, t):
    """
    Find all distinct k-mers forming (L, t)-clumps in Genome.

    :param text: DNA
    :param k: length of kmers
    :param L: ori length
    :param t: minimum number of occurences of kmer in L
    :return:
    """
    all_kmers = set()
    for start in range(0, len(text) - L + 1):
        kmers, max_cnt = frequent_kmers(text[start:start + L], k)
        if max_cnt >= t:
            all_kmers.update(kmers)

    return all_kmers


def get_keys_(dict_, t):
    return {k for (k,v) in dict_.items() if v >= t}

def find_clumps_opti(text, k, L, t):
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
