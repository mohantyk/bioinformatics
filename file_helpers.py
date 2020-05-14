from collections import defaultdict

## File Read Helpers
def read_list(list_with_spaces, dtype=int):
    return [dtype(x) for x in list_with_spaces.split()]

def read_matrix(data):
    matrix = []
    for line in data.splitlines():
        matrix.append(read_list(line))
    return matrix

def get_data(filename):
    with open(filename) as f:
        data = f.readlines()
    return data

def read_adjacency(lines):
    adjacency = {}
    for line in lines:
        line = line.strip()
        (head, tails) = line.split(' -> ')
        tails = tails.split(',')
        adjacency[int(head)] = [int(node) for node in tails]
    return adjacency

def read_weighted_adjacency(graph):
    adjacency = defaultdict(list)
    for line in graph.splitlines():
        node, nxt = line.split('->')
        nxt_node, weight = nxt.split(':')
        adjacency[int(node)].append((int(nxt_node), int(weight)))
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

#-------------------------------------------------------------------------------
# Converters

def convert_txt_to_csv(filename):
    '''
    Converts a space separated columnar text file to a csv file
    '''
    csv_file = filename.replace('.txt', '.csv')
    csv_lines = []
    with open(filename, 'r') as f:
        for line in f:
            csv_lines.append(','.join(line.split()))
    with open(csv_file, 'w') as f:
        f.write('\n'.join(csv_lines))


def char2int(char):
    return ord(char) - ord('a')

def int2char(num):
    return chr(97+num)

### File Write Helpers

def adjacency_to_file(adjacency, filename):
    with open(filename, 'w') as f:
        for k, v in adjacency.items():
            f.write(f'{k} -> {", ".join(v)}\n')

def weighted_adjacency_to_file(adjacency, filename):
    nodes = sorted(adjacency.keys())
    with open(filename, 'w') as f:
        for node in nodes:
            for nghbr, wght in adjacency[node].items():
                print(f'{node}->{nghbr}:{int(wght)}', file=f)

def path_to_file(path, filename):
    path_str = '->'.join(str(node) for node in path)
    with open(filename, 'w+') as f:
        f.write(path_str+'\n')
