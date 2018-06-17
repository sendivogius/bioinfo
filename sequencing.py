import itertools

import motifs
from printing import print_list_as_path


def composition(dna, k):
    return sorted(motifs.get_all_kmers(dna, k, ordered=True))


def genome_from_sequence(kmers_sequence):
    k = len(kmers_sequence[0])
    return ''.join([l[0] for l in kmers_sequence[:(len(kmers_sequence)-1)]])+kmers_sequence[-1]


def print_formatted(d):
    for k in sorted(d.keys()):
        print(k, '->', ','.join(sorted(d[k])))


def graph_from_sequence(kmers_sequence):
    k = len(kmers_sequence[0])
    sorted_kmers = sorted(kmers_sequence)
    all_kmers = {kmer1: [kmer2 for kmer2 in sorted_kmers if kmer1[1:k] == kmer2[:k - 1]] for kmer1 in sorted_kmers}
    v = {k: v for k, v in all_kmers.items() if v}
    return v

def kmer_prefix(dnas):
    if type(dnas) is str:
        return dnas[:-1]
    return tuple((d[:-1] for d in dnas))

def kmer_suffix(dnas):
    if type(dnas) is str:
        return dnas[1:]
    return tuple((d[1:] for d in dnas))

def debrujin_from_kmers(kmers):
    graph = dict()
    for kmer1, kmer2 in ((kmer_prefix(kmer), kmer_suffix(kmer)) for kmer in kmers):
        if kmer1 in graph:
            graph[kmer1].append(kmer2)
        else:
            graph[kmer1] = [kmer2]
    for key, val in graph.items():
        graph[key] = sorted(val)
    # print_formatted(graph)
    return graph


def debrujin_from_dna(dna, k):
    kmers = motifs.get_all_kmers(dna, k, ordered=True)
    return debrujin_from_kmers(kmers)


def eulerian_cycle(graph):
    nodes_stack = [next(iter(graph))]
    euler_cycle = []
    while nodes_stack:
        top_node = nodes_stack[-1]
        if graph[top_node]:
            nodes_stack.append(graph[top_node].pop())
        else:
            euler_cycle.append(nodes_stack.pop())
    euler_cycle.reverse()
    return euler_cycle


def _get_terminal_nodes(graph):
    out_degrees = {k: len(v) for k, v in graph.items()}
    all_edges = [el for edges in graph.values() for el in edges]
    in_degrees = { k: all_edges.count(k) for k in set(all_edges)}
    degrees = {key: (out_degrees[key] - in_degrees.get(key, 0)) for key in out_degrees.keys()}
    out_node = [k for k, v in degrees.items() if v == -1]
    if out_node:
        out_node = out_node[0]
    else:
        out_node = next(iter(set(in_degrees.keys()) - set(out_degrees.keys())))
    in_node = [k for k, v in degrees.items() if v == 1]
    return in_node[0], out_node


def eulerian_path(graph):
    start_node, end_node = _get_terminal_nodes(graph)
    graph[end_node] = graph.get(end_node, []) + [start_node]
    cycle = eulerian_cycle(graph)[1:]
    #breaking the cycle at added dummy edge
    break_position = next(i for i in range(len(cycle)-1) if cycle[i:i+2] == [end_node, start_node])+1
    path = cycle[break_position:] + cycle[:break_position]
    # print_list_as_path(path)
    return path


def string_reconstruction(kmers):
    graph = debrujin_from_kmers(kmers)
    path = eulerian_path(graph)
    return genome_from_sequence(path)


def universal_string(k):
    kmers = list(''.join(t) for t in itertools.product('01', repeat=k))
    graph = debrujin_from_kmers(kmers)
    path = eulerian_cycle(graph)
    return genome_from_sequence(path)[:-k+1]

def paired_composition(dna, k, d):
    pairs = [(dna[i:i+k], dna[i+k+d:i+2*k+d]) for i in range(len(dna)-2*k-d+1)]
    return sorted(pairs)

