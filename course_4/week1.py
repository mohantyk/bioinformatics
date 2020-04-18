import numpy as np

def distances(n, adjacency):
    '''
    inputs:
        n: n is the number of leaf nodes (0...n-1)
        adjacency : weighted adjacency matrix node: [(nghbr0, weight0), (nghbr1, weight1)]
    outputs:
        distance matrix
    '''
    # Create initial distance graph - single edges
    max_node = max(adjacency.keys())
    distances = np.zeros((max_node+1, max_node+1), int)
    for node in adjacency:
        for (nghbr, weight) in adjacency[node]:
            distances[node, nghbr] = weight

    # Iteratively remove node and update weights to its neighbors
    for node in reversed(range(n, max_node+1)):
        nghbrs = np.where(distances[node,:])[0] # Find all its neighbors
        for nghbr in nghbrs:
            for idx in range(node):
                if distances[nghbr, idx] == 0 and idx in nghbrs: # Update only if weight is unknown
                    distances[nghbr, idx] = distances[node, idx] + distances[nghbr, node]
        distances = distances[:node, :node] # Remove the node
        np.fill_diagonal(distances, 0)

    return distances
