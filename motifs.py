import functools
import itertools
import operator
import collections
from math import log2
from random import randrange
import numpy.random

from utils import is_dna

_rev_mapping = {'C': 'G', 'G': 'C', 'A': 'T', 'T': 'A'}


def reverse_complement_strand(dna):
    """Compute reverse complement strand for given dna (1C in BA_AALA)
    :param dna: ACTG string
    :return: reversed complement strand
    """
    assert (is_dna(dna))
    return ''.join(_rev_mapping[nn] for nn in dna[::-1])


def get_all_kmers(pattern, k, ordered=False):
    """Returns all kmers for pattern
    Depending on `ordered` returns all kmers as they appear in `pattern` or set of unique kmers

    :param pattern: string to look for kmers
    :param k: length of kmer
    :param ordered: boolean flag
    :return: list or se of kmers, depending o `ordered` parameter
    """
    ordered_kmers = [pattern[i:i + k] for i in range(len(pattern) - k + 1)]
    if ordered:
        return ordered_kmers
    return set(ordered_kmers)


def find_occurrences(text, pattern, d=0):
    """Get indexes of (possibly overlapping) occurrences of `pattern` in `text` with up to `d` mismatches
    1D (k=0) and 1H (k>0) in BA_AALA
    :param text: string to search in
    :param pattern: string to look for
    :param d: int, maximum number of mismatches (aka Hamming distance)
    :return: list of indices where pattern appears in text
    """
    idx_of_last_pattern = len(text) - len(pattern)
    return [i for i in range(idx_of_last_pattern + 1) if hamming(text[i:i + len(pattern)], pattern) <= d]


def count_occurrences(text, pattern, d=0):
    """Count number of (possibly overlapping) occurrences of pattern in text (1A in BA_AALA)
    :param text: string to search in
    :param pattern: string to look for
    :param d: int, maximum number of mismatches (aka Hamming distance)
    :return: int, number of times `pattern` appears in `text`
    """
    return len(find_occurrences(text, pattern, d))


def frequent_kmers(genome, k, d=0, reverse=False):
    """Get most frequent kmers in string with up to `d` mismatches and reverse complement
    Problem 1B (d=0, reverse=False), 1I (d>1, reverse=True) and 1J (revere=True) in BA_AALA

    :param genome: ACTG string
    :param k: length of kmer
    :param d: number of allowed mismatches
    :param reverse: flag whether include reverse complement strand
    :return: tuple ({kmer1,kmer2,..), num_of_occurences
    """
    assert (is_dna(genome))

    max_kmers = []
    max_kmers_cnt = 0

    reversed = reverse_complement_strand(genome) if reverse else ''

    kmers = get_all_kmers(genome, k)
    kmers_neighbourhood = {p for kmer in kmers for p in get_neighbours(kmer, d)}
    for i, kmer in enumerate(kmers_neighbourhood):
        kmer_counts = count_occurrences(genome, kmer, d) + count_occurrences(reversed, kmer, d)
        if kmer_counts == max_kmers_cnt:
            max_kmers.append(kmer)
        elif kmer_counts > max_kmers_cnt:
            max_kmers = [kmer]
            max_kmers_cnt = kmer_counts
        # if i % 5000 == 0:
        #     print(i)
    return set(max_kmers), max_kmers_cnt


def _get_neighbours(kmer):
    """Get immediate neighbours of kmer

    :param kmer: ACTG string
    :return: set of neighbours kmers
    """
    assert (is_dna(kmer))
    bases = 'ACTG'
    result = set()
    for i in range(len(kmer)):
        for base in bases:
            result.add(kmer[:i] + base + kmer[(i + 1):])
    return result


def get_neighbours(kmer, max_d):
    """Generate all kmers of distance at most max_d to kmer (1N in BA_AALA)

    :param kmer: ACTG string
    :param max_d: maximum Hamming distance
    :return: set of neighbours dnas
    """
    assert (is_dna(kmer))
    result = set([kmer])
    for i in range(max_d):
        addded = set()
        for kmer in result:
            addded |= _get_neighbours(kmer)
        result |= addded
    return result


def _get_keys(dict_, t):
    """Gets keys for items having value at least t"""
    return {k for (k, v) in dict_.items() if v >= t}


def find_clumps(genome, k, L, t):
    """Find all distinct k-mers forming (L, t)-clumps in genome (Problem 1E in BA_AALA)
    (L, t)-clump is an interval of length L in which k-mer appears at least t times
    :param genome: ACTG string
    :param k: length of kmer
    :param L: ori length
    :param t: minimum number of occurrences of kmer in L
    :return: set of distinct kmers forming (L-T)-clumps
    """
    assert (is_dna(genome))
    counts = collections.defaultdict(int)

    # compute counts of kmers in first L-length part of genome
    for k_start in range(L - k + 1):
        counts[genome[k_start:k_start + k]] += 1
    kmers = _get_keys(counts, t)

    # slide L-length window and update counts
    # remove previous leftmost kmer and add new kmer being rightmost in current window
    for L_start in range(1, len(genome) - L + 1):
        counts[genome[L_start - 1:L_start + k - 1]] -= 1
        new_kmer = genome[L_start + L - k:L_start + L]
        counts[new_kmer] += 1
        if counts[new_kmer] >= t:
            kmers.add(new_kmer)
    return kmers


def pattern_to_number(dna):
    """Computes number corresponding to dna pattern (1L in BA_AALA)

    :param dna: ACTG string
    :return: int
    """
    assert (is_dna(dna))
    idx = 'ACGT'
    return sum(idx.index(dna_base) * 4 ** i for i, dna_base in enumerate(dna[::-1]))


def number_to_pattern(number, length):
    """Computes dna pattern for given number (1K in BA_AALA)

    :param number:
    :param length: length of returned kmer
    :return: kmer for which pattern_to_number returns `number`
    """

    idx = 'ACGT'
    pattern = ''
    while number > 0:
        pattern += idx[number % 4]
        number //= 4
    return idx[0] * (length - len(pattern)) + pattern[::-1]


def get_skew(genome):
    """Computes skew of genome
    Skew is defined as running total of differences between G and C
    :param genome: ACTG string
    :return: list of same length as `genome` with skew values
    """
    assert (is_dna(genome))
    skew = [0] * (len(genome) + 1)
    for i, base in enumerate(genome, 1):
        if base == 'C':
            skew[i] = skew[i - 1] - 1
        elif base == 'G':
            skew[i] = skew[i - 1] + 1
        else:
            skew[i] = skew[i - 1]
    return skew


def get_min_skew_position(genome):
    """Find a position in a genome where the skew attain a minumum (1F in BA_AALA)

    :param genome: ACTG string
    :return: list of all integers minimizing skew
    """
    assert (is_dna(genome))
    skew = get_skew(genome)
    min_skew = min(skew)
    return [pos for (pos, sk) in enumerate(skew) if sk == min_skew]


def plot_skew(genome):
    assert (is_dna(genome))
    import matplotlib.pyplot as plt
    skew = get_skew(genome)
    plt.plot(skew, '-')
    plt.title('Skew plot')
    plt.ylabel('skew')
    plt.xlabel('position')
    plt.show()


def hamming(text1, text2):
    """Calculate Hamming distance between two strings of equal length (1G in BA_AALA)

    :return: The Hamming distance between input strings
    """
    assert (len(text1) == len(text2))
    return sum(a != b for a, b in zip(text1, text2))


def motif_enumeration(dnas, k, d):
    """Brute force algorithm for solving `Implanted Motif Problem` (2A in BA_AALA)
    (k,d)-motif is a string is a kmer that appears in every dna string with at most d mismatches

    :param dnas: list of ACTG strings
    :param k: length of kmer
    :param d: number of mismatches
    :return: set of (k,d)-motifs
    """
    assert(all(is_dna(dna) for dna in dnas))
    patterns = set()
    for kmer in get_all_kmers(dnas[0], k):
        for kmer_neigh in get_neighbours(kmer, d):
            counts = [count_occurrences(dna, kmer_neigh, d) for dna in dnas]
            if all(counts):
                patterns.add(kmer_neigh)
    return patterns


def get_profile(motifs, relative=True, pseudocounts=False):
    """Computes profile of motifs array.

    :param motifs: list of ACTG strings
    :param relative: flag for computing frequencies instead of counts
    :param pseudocounts: flag for using Laplacian smoothing
    :return: list of dicts with counts or frequencies
    """
    assert (all(is_dna(motif) for motif in motifs))
    zeroes = dict(itertools.zip_longest('ACTG', [], fillvalue=0))
    counts = [{**zeroes, **dict(collections.Counter(nucleotides_at_pos))} for nucleotides_at_pos in zip(*motifs)]

    if pseudocounts:
        counts = [{nuc: cnt + 1 for nuc, cnt in pos_cnts.items()} for pos_cnts in counts]

    if relative:
        return [{pos: cnt / sum(pos_cnts.values()) for pos, cnt in pos_cnts.items()} for pos_cnts in counts]

    return counts


def get_consensus_string(motifs):
    """Get most probable kmer describes by motifs

    :param motifs: list of ACTG strings
    :return:
    """
    profile = get_profile(motifs)
    return ''.join(max(position.items(), key=operator.itemgetter(1))[0] for position in profile)


def score_motifs(motifs, entropy=False):
    """Score motifs

    :param motifs: list of ACTG strings
    :param entropy: flag whether use entropy or number of mismatches
    :return: cumulative entropy or number of mistmatches
    """
    assert (all(is_dna(motif) for motif in motifs))

    def _calc_entropy(values):
        return sum(-v * log2(v) if v else 0.0 for v in values)

    if entropy:
        profile = get_profile(motifs, relative=True)
        return sum(_calc_entropy(val.values()) for val in profile)

    profile = get_profile(motifs, relative=False)
    consensus = get_consensus_string(motifs)
    return sum((sum({**p, **{c: 0}}.values()) for c, p in zip(consensus, profile)))


def dist(pattern, dnas):
    """Computes distance between pattern and dnas as sum of distances of pattern to every dna
    Distance between pattern and dna is minimum Hamming distance between pattern and every kmer in dna.

    :param pattern: ACTG string
    :param dnas: ACTG string or list of ACTG strings
    :return:
    """
    if type(dnas) == str:
        dnas = [dnas]

    assert(is_dna(pattern))
    assert(all(is_dna(dna) for dna in dnas))

    def dist_to_single_dna(pat, dna_string):
        return min(hamming(pat, kmer) for kmer in get_all_kmers(dna_string, len(pat)))

    return sum(dist_to_single_dna(pattern, dna) for dna in dnas)


def median_string(dnas, k):
    """Find a median string (2B in BA_AALA)
    Median string is a kmer that minimizes d(pattern, dnas) among all possible kmers

    :param dnas: list of ACTG strings
    :param k: length of kmer
    :return: tuple(median_string, dist_median_dnas)
    """
    assert(all(is_dna(dna) for dna in dnas))
    all_kmers = map(''.join, itertools.product('ACTG', repeat=k))
    return min(((kmer, dist(kmer, dnas)) for kmer in all_kmers), key=operator.itemgetter(1))


def score_kmer_on_profile(kmer, profile):
    """Calculates probability that profile generates kmer

    :param kmer: ACTG string
    :param profile: list of dicts with probabilities
    :return:
    """
    assert(is_dna(kmer))
    assert (len(kmer) == len(profile))
    proba_kmer = (p[k] for p, k in zip(profile, kmer))
    return functools.reduce(operator.mul, proba_kmer, 1)


def most_probable_kmer_from_profile(dna, profile):
    """Find a most probably kmer in a string (2C in BA_AALA)

    :param dna: ACTG string
    :param profile: list of dicts with probabilities
    :return: profile most probably kmer
    """
    assert(is_dna(dna))

    k = len(profile)
    score, kmer = max([(score_kmer_on_profile(kmer, profile), kmer) for kmer in get_all_kmers(dna, k, ordered=True)],
               key=operator.itemgetter(0))
    return kmer


def greedy_motif_search(dnas, k):
    """Finds bets motifs in greedy way (selecting best motifs at each iteration)
    Includes pseudocounts (Laplacian smoothing)
    (2D and 2E in BA_AALA)

    :param dnas: list of ACTG strings
    :param k: length of kmer
    :return: tuple(best_motifs, best_score)
    """
    assert (all(is_dna(dna) for dna in dnas))

    best_motifs, best_score = _get_first_kmers(dnas, k)

    for kmer in get_all_kmers(dnas[0], k, ordered=True):
        motifs = [kmer]
        for i in range(1, len(dnas)):
            profile = get_profile(motifs, pseudocounts=True)
            motif = most_probable_kmer_from_profile(dnas[i], profile)
            motifs.append(motif)
        motifs_score = score_motifs(motifs)
        if motifs_score < best_score:
            best_motifs = motifs
            best_score = motifs_score
    return best_motifs, best_score


def _get_first_kmers(dnas, k):
    """Get first kmer from string as motifs

    :param dnas: list of ACTG strings
    :param k: length of kmer
    :return: (motifs, motifs_score)
    """

    best_motifs = [dna[:k] for dna in dnas]
    best_score = score_motifs(best_motifs)
    return best_motifs, best_score


def get_random_motif(dna, k):
    """Get random kmer from dna

    :param dna: ACTG string
    :param k: length of kmer
    :return: random kmer from dns
    """
    assert(is_dna(dna))
    assert(len(dna) - k + 1 > 0)
    start = randrange(len(dna) - k + 1)
    return dna[start:start + k]


def randomized_motif_search(dnas, k, times=2000):
    """Finds bets motifs by randomly selecting set of initial kmers and then iteratively improving profile
    (2F in BA_AALA)

    :param dnas: list of ACTG strings
    :param k: length of kmer
    :param times: number of time the algorithm will be run with different initial kmers
    :return: tuple(best_motifs, best_score) from `times` runs
    """
    assert (all(is_dna(dna) for dna in dnas))

    best_motifs, best_score = _get_first_kmers(dnas, k)
    for run in range(times):
        motifs = [get_random_motif(d, k) for d in dnas]
        while True:
            profile = get_profile(motifs, pseudocounts=True)
            motifs = [most_probable_kmer_from_profile(dna, profile) for dna in dnas]
            motifs_score = score_motifs(motifs)
            if motifs_score < best_score:
                best_motifs = motifs
                best_score = motifs_score
            else:
                break
    return best_motifs, best_score


def get_profile_randomly_generated_kmer(dna, profile):
    """Get random kmer from dna with chance of being selected defined by profile

    :param dna: ACTG string
    :param profile: list of dicts with probabilities
    :return: random kmer with probabilities defined by profile
    """
    assert(is_dna(dna))

    k = len(profile)
    kmer_probs = [(kmer, score_kmer_on_profile(kmer, profile)) for kmer in get_all_kmers(dna, k)]
    kmers, probs = zip(*kmer_probs)
    total_prob = sum(probs)
    normalized_probs = [p / total_prob for p in probs]
    draw = numpy.random.choice(kmers, 1, p=normalized_probs)
    return draw[0]


def gibbs_sampler(dnas, k, N=100, times=1000):
    """Finds bets motifs by using Gibbs sampling (2G in BA_AALA)
    Excude single kmer from current motifs set and decide whether to keep it or replace.


    :param dnas: list of ACTG strings
    :param k: length of kmer
    :param N: number of iterations in one sampling run
    :param times: number of time sampling will be performed with different initial kmers
    :return: tuple(best_motifs, best_score) from `times` runs
    """
    assert (all(is_dna(dna) for dna in dnas))

    best_motifs, best_score = _get_first_kmers(dnas, k)
    n_dna = len(dnas)
    for run in range(times):
        motifs = [get_random_motif(dna, k) for dna in dnas]
        for _ in range(N):
            exclude = randrange(n_dna)
            motifs_excluded = [motifs[i] for i in range(n_dna) if i != exclude]
            profile = get_profile(motifs_excluded, pseudocounts=True)
            motifs[exclude] = get_profile_randomly_generated_kmer(dnas[exclude], profile)
            motifs_score = score_motifs(motifs)
            if motifs_score < best_score:
                best_motifs = motifs
                best_score = motifs_score
    return best_motifs, best_score


if __name__ == "__main__":
    ecoli = open('data\\E_coli.txt').readline()
    plot_skew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')
