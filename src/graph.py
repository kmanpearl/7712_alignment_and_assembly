import datetime

import pandas as pd


# find edges between nodes
def create_graph(kmers):
    edges = []
    for k1 in range(len(kmers)):
        for k2 in range(1, len(kmers)):
            if kmers[k1].sufix == kmers[k2].prefix:
                edges.append([kmers[k1], kmers[k2]])
    print(f"{datetime.datetime.now()}: found all edges of the graph")
    return edges


# turn edges into an adjaceny matrix
def create_adjacency_matrix(kmers, edges, save):
    nodes = []
    for kmer in kmers:
        if kmer.sequence not in nodes:
            nodes.append(kmer.sequence)
    adj_matrix = pd.DataFrame(0, index=nodes, columns=nodes)
    for edge1, edge2 in edges:
        source = edge1.sequence
        target = edge2.sequence
        adj_matrix.loc[source, target] = 1
    adj_matrix = adj_matrix.astype(int)
    print(f"{datetime.datetime.now()}: created adjacency matrix")
    if save:
        adj_matrix.to_csv("output/adjacency_matrix.csv")
    return adj_matrix


def find_start_stop_nodes(adj_matrix):
    """
    Finds all possible start and stop nodes in a directed graph represented by an adjacency matrix.

    Parameters:
        adj_matrix (pandas.DataFrame): The adjacency matrix of the graph.

    Returns:
        A tuple containing two lists: a list of all possible start nodes, and a list of all possible
        stop nodes. Each list contains the node labels as strings.
    """
    # initialize sets of start and stop nodes
    start_nodes = list(adj_matrix.columns)
    stop_nodes = list(adj_matrix.columns)

    # iterate over the rows and columns of the adjacency matrix
    for node in adj_matrix.index:
        if any(adj_matrix.loc[:, node] != 0):
            start_nodes.remove(node)
        if any(adj_matrix.loc[node, :] != 0):
            stop_nodes.remove(node)

    # convert the sets to lists and return the result
    print(f"{datetime.datetime.now()}: found all start and stop nodes")
    return list(start_nodes), list(stop_nodes)


class Paths:
    def __init__(self, contig_id, start_node, stop_node, path):
        self.contig_id = contig_id
        self.start_node = start_node
        self.stop_node = stop_node
        self.path = path


def find_all_paths(adj_matrix):
    """
    Finds all possible paths in a directed graph represented by an adjacency matrix.

    Parameters:
        adj_matrix (pandas.DataFrame): The adjacency matrix of the graph.

    Returns:
        A dictionary of all possible paths, keyed by the tuple (start_node, stop_node), where
        start_node and stop_node are the labels of the start and stop nodes, respectively. Each value
        in the dictionary is a list of lists, where each inner list represents a path from the start
        node to the stop node.
    """
    # find all possible start and stop nodes
    start_nodes, stop_nodes = find_start_stop_nodes(adj_matrix)

    # initialize the result dictionary
    all_paths = []
    contig_id = 0
    # iterate over all pairs of start and stop nodes
    for start_node in start_nodes:
        for stop_node in stop_nodes:
            # skip pairs where the start and stop nodes are the same
            if start_node == stop_node:
                continue

            visited = set()
            stack = [(start_node, [start_node])]
            paths = []

            while stack:
                (node, path) = stack.pop()
                if node == stop_node:
                    paths.append(path)
                    continue
                if node in visited:
                    continue
                visited.add(node)
                for neighbor in adj_matrix.loc[node, :][
                    adj_matrix.loc[node, :] != 0
                ].index:
                    stack.append((neighbor, path + [neighbor]))

            # add the paths to the result dictionary
            if paths:
                for path in paths:
                    p = Paths(
                        contig_id=contig_id,
                        start_node=start_node,
                        stop_node=stop_node,
                        path=path,
                    )
                    all_paths.append(p)
                    contig_id += 1
    print(f"{datetime.datetime.now()}: found {len(all_paths)} paths through the graph")
    return all_paths


def graph_traversal(kmers, save):
    edges = create_graph(kmers)
    adj_matrix = create_adjacency_matrix(kmers, edges, save)
    all_paths = find_all_paths(adj_matrix)
    return all_paths
