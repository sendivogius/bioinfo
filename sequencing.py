import motifs


def composition(dna, k):
    return sorted(motifs.get_all_kmers(dna, k, ordered=True))


def genome_from_sequence(kmers_sequence):
    k = len(kmers_sequence[0])
    return (''.join(kmers_sequence[::k])) + kmers_sequence[-1][len(kmers_sequence) // k:]


def print_formatted(d):
    for k in sorted(d.keys()):
        print(k, '->', ','.join(d[k]))


def graph_from_sequence(kmers_sequence):
    k = len(kmers_sequence[0])
    all_kmers = {kmer1: {kmer2 for kmer2 in kmers_sequence if kmer1[1:k] == kmer2[:k - 1]} for kmer1 in kmers_sequence}
    v = {k: v for k, v in all_kmers.items() if v}
    print_formatted(v)
    # return v
