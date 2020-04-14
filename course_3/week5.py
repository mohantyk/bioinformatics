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
