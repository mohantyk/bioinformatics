def greedy_sort_permutation(perm):
    n = len(perm)
    identity = list(range(1, n+1))

    last = perm
    steps = []
    while last != identity:
        for i in range(1, n+1):
            if last[i-1] != i:
                break
        k = i
        for i in range(k, n+1):
            if last[i-1] in (k, -k):
                break
        l = i
        reversal = list(reversed([-x for x in last[k-1:l]]))
        new = last[:k-1] + reversal + last[l:]
        steps.append(new)
        last = new

    return steps


def num_breakpoints(perm):
    count = 0
    n = len(perm)
    extended = [0] + perm + [n+1]
    for i in range(n+1):
        k = extended[i]
        l = extended[i+1]
        if l-k != 1:
            count += 1
    return count
