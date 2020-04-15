def chromosome_to_cycle(perm):
    nodes = []
    for n in perm:
        head = 2*abs(n)
        tail = 2*abs(n)-1
        if n > 0:
            nodes.append(tail)
            nodes.append(head)
        else:
            nodes.append(head)
            nodes.append(tail)
    return nodes

def pairwise(iterable):
    a = iter(iterable)
    return zip(a, a)

def cycle_to_chromosome(cycle):
    perm = []
    for x, y in pairwise(cycle):
        if y > x:
            perm.append(y//2)
        else:
            perm.append(-x//2)
    return perm

def colored_edges(genome):
    edges = []
    for chromosome in genome:
        nodes = chromosome_to_cycle(chromosome)
        for edge in pairwise(nodes[1:] + [nodes[0]]):
            edges.append(edge)
    return edges

def cycles_in_genome_graph(edges):
    cycles = []
    end_node = None
    for edge in edges:
        if end_node is None:
            cycle = []
            start_node = edge[0]
            if start_node%2 == 0:
                end_node = start_node - 1
            else:
                end_node = start_node + 1
        cycle.append(edge)
        if edge[-1] == end_node:
            end_node = None
            cycles.append(cycle)
    return cycles


def graph_to_genome(edges):
    raise NotImplementedError
