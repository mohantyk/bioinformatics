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
