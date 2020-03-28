
## File Read Helpers

def read_adjacency( filename ):
    with open(filename, 'r') as f:
        lines = f.readlines()
    adjacency = {}
    for line in lines:
        line = line[:-1]
        (head, tails) = line.split(' -> ')
        tails = tails.split(',')
        adjacency[int(head)] = [int(node) for node in tails]
    return adjacency


def transform_pairs_to_tuples(paired_read):
    return tuple(paired_read.split('|'))


def read_pairs(filename):
    with open(filename, 'r') as f:
        data = f.readlines()
    k,d = data[0].strip().split()
    k = int(k)
    d = int(d)
    kmer_pairs = [transform_pairs_to_tuples(line.strip()) for line in data[1:]]
    return k, d, kmer_pairs

### File Write Helpers

def adjacency_to_file(adjacency, filename):
    with open(filename, 'w') as f:
        for k, v in adjacency.items():
            f.write(f'{k} -> {", ".join(v)}\n')

def path_to_file(path, filename):
    path_str = '->'.join(str(node) for node in path)
    with open(filename, 'w+') as f:
        f.write(path_str+'\n')
