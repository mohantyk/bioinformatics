from collections import deque

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
    genome = []
    for cycle in cycles_in_genome_graph(edges):
        nodes = deque(node for edge in cycle for node in edge)
        nodes.rotate()
        chromosome = cycle_to_chromosome(nodes)
        genome.append(chromosome)
    return genome

def two_break_on_genome_graph(genome_graph, i1, i2, i3, i4):
    '''
    Convert (i1, i2) -> (i4, i2)
    and     (i3, i4) -> (i3, i1)
    '''
    new_graph = []
    for edge in genome_graph:
        new_edge = edge
        if edge == (i1, i2):
            new_edge = (i4, i2)
        elif edge == (i2, i1):
            new_edge = (i2, i4)
        elif edge == (i3, i4):
            new_edge = (i3, i1)
        elif edge == (i4, i3):
            new_edge = (i1, i3)
        new_graph.append(new_edge)
    return new_graph
