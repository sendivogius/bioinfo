import itertools

from baal import motifs


def composition(dna, k):
    return sorted(motifs.get_all_kmers(dna, k, ordered=True))


def genome_from_sequence(kmers_sequence):
    return ''.join([l[0] for l in kmers_sequence[:(len(kmers_sequence) - 1)]]) + kmers_sequence[-1]


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
    return tuple(d[:-1] for d in dnas)


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
    out_degrees = {k: _get_out_degree(graph, k) for k in graph}
    in_degrees = {k: _get_in_degree(graph, k) for k in graph}
    degrees = {key: (out_degrees[key] - in_degrees.get(key, 1)) for key in out_degrees.keys()}
    out_node = [k for k, v in degrees.items() if v == -1]
    if out_node:
        out_node = out_node[0]
    else:
        all_edges = {el for edges in graph.values() for el in edges}
        out_node = next(iter(all_edges - set(out_degrees.keys())))
    in_node = [k for k, v in degrees.items() if v == 1]
    return in_node[0], out_node


def eulerian_path(graph):
    start_node, end_node = _get_terminal_nodes(graph)
    # add dummy edge
    graph[end_node] = graph.get(end_node, []) + [start_node]
    cycle = eulerian_cycle(graph)[1:]
    # breaking the cycle at added dummy edge
    break_position = next(i for i in range(len(cycle) - 1) if cycle[i:i + 2] == [end_node, start_node]) + 1
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
    return genome_from_sequence(path)[:-k + 1]


def paired_composition(dna, k, d):
    pairs = [(dna[i:i + k], dna[i + k + d:i + 2 * k + d]) for i in range(len(dna) - 2 * k - d + 1)]
    return sorted(pairs)


def genome_from_pair_sequence(seq, d):
    def match(gene, pattern):
        return all((g == '_' or g == p for (g, p) in zip(gene, pattern)))

    k = len(seq[0][0])
    genome = ['_'] * (len(seq) + 2 * k + d - 1)
    for i, (k1, k2) in enumerate(seq):
        k1_start = i
        k2_start = k1_start + k + d
        # check if k1 and k2 can be added to already constructed genome
        if match(genome[k1_start:k1_start + k], k1) and match(genome[k2_start:k2_start + k], k2):
            genome[k1_start:k1_start + k] = k1
            genome[k2_start:k2_start + k] = k2
        else:
            return None

    return ''.join(genome)


def string_pair_reconstruction(kmers, d):
    graph = debrujin_from_kmers(kmers)
    path = eulerian_path(graph)
    genome = genome_from_pair_sequence(path, d + 1)
    return genome
    # for path in eulerian_path(graph):
    #     reconstructed_path = genome_from_pair_sequence(path, d + 1)
    #     if reconstructed_path:
    #         return reconstructed_path
    # return None


def maximal_nonbranching_paths(graph):
    def one_in_one_out(graph, node):
        out_degree = _get_out_degree(graph, node)
        in_degree = _get_in_degree(graph, node)
        return (1 == in_degree == out_degree)

    starting_nodes = (node for node in graph if not one_in_one_out(graph, node) and _get_out_degree(graph, node) > 0);
    for node in starting_nodes:
        # try each outgoing edges
        for w in graph[node]:
            non_branching_path = [node, w]
            while one_in_one_out(graph, w):
                w = graph[w][0]
                non_branching_path.append(w)
            yield non_branching_path

    possible_cycle_nodes = (node for node in graph if one_in_one_out(graph, node))
    for node in possible_cycle_nodes:
        cycle = [node]
        next_node = graph[node][0]
        while one_in_one_out(graph, next_node):
            cycle.append(next_node)
            next_node = graph[next_node][0]
            if cycle[0] == cycle[-1]:
                yield cycle
                for n in cycle:
                    graph[n] = []
                break


def _get_in_degree(graph, node):
    all_edges = [el for edges in graph.values() for el in edges]
    return all_edges.count(node)


def _get_out_degree(graph, node):
    return len(graph.get(node, []))


def get_contigs(kmers):
    graph = debrujin_from_kmers(kmers)
    mnbp = maximal_nonbranching_paths(graph)
    return [genome_from_sequence(s) for s in mnbp]
