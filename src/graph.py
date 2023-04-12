__author__ = "Keenan Manpearl"
__date__ = "2023/04/09"

"""
Functions needed to assemble de bruijn graphs   
"""

import datetime

import pandas as pd


class Paths:
    """
    holds information about the path through the graph between a start and stop node
    """

    def __init__(self, contig_id, start_node, stop_node, path):
        self.contig_id = contig_id
        self.start_node = start_node
        self.stop_node = stop_node
        self.path = path


# find edges between nodes
def create_graph(kmers):
    """
    create a de bruijin graph

    Args:
        kmers (list): instances of class Kmer that make up all reads

    Returns:
        list: list of touples containing a source node and a target node
    """
    edges = []
    for k1 in range(len(kmers)):
        for k2 in range(1, len(kmers)):
            if kmers[k1].sufix == kmers[k2].prefix:
                edges.append([kmers[k1], kmers[k2]])
    print(f"{datetime.datetime.now()}: found all edges of the graph")
    return edges


# turn edges into an adjaceny matrix
def create_adjacency_matrix(kmers, edges, save, out_dir):
    """
    turn edge list into an adjacency matrix

    Args:
        kmers (list): instances of class Kmers present in the reads
        edges (list): touples of start and stop nodes
        save (bool): whether to save matrix as a csv

    Returns:
        pandas.DataFrame: adjacency matrix where source node are rows and target nodes are columns
    """
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
        adj_matrix.to_csv(f"{out_dir}/adjacency_matrix.csv")
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
    if len(start_nodes) == 0:
        raise Exception("No start nodes found, graph is cyclic")
    if len(stop_nodes) == 0:
        raise Exception("No stop nodes found, graph is cyclic")

    # convert the sets to lists and return the result
    print(f"{datetime.datetime.now()}: found all start and stop nodes")
    return list(start_nodes), list(stop_nodes)


def find_all_paths(adj_matrix):
    """
    Finds all possible paths in a directed graph.

    Parameters:
        adj_matrix (pandas.DataFrame): The adjacency matrix of the graph.

    Returns:
        list : instance of class Path representing all possible paths between start and stop nodes
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
                    # if we get to a cycle, discard the path
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


def graph_traversal(kmers, save, out_dir):
    """
    wrapper function for creating and traversing graph

    Args:
        kmers (list): instances of class Kmer represening all read kmers
        save (bool): whether to save adjacency matrix as csv file

    Returns:
        list: instance of class Path that stores all possible paths through the graph
    """
    edges = create_graph(kmers)
    adj_matrix = create_adjacency_matrix(kmers, edges, save, out_dir)
    all_paths = find_all_paths(adj_matrix)
    return all_paths
