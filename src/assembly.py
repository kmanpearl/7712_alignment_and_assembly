__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

"""
Functions needed to assemble sequence reads 
"""

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


def get_contig_kmers(all_paths, read_kmers):
    contig_kmers = {}
    for path in all_paths:
        contig_id = path.contig_id
        for node in path.path:
            for kmer in read_kmers:
                if node == kmer.sequence:
                    if contig_id in contig_kmers.keys():
                        contig_kmers[contig_id].append(kmer)
                    else:
                        contig_kmers[contig_id] = [kmer]
                    break
    return contig_kmers


class Contig:
    def __init__(self, contig_id, aligned_reads, sequence, direction):
        self.contig_id = contig_id
        self.sequence = sequence
        self.direction = direction
        self.aligned_reads = aligned_reads

    def add_kmers(self, kmers):
        self.kmers = kmers


class AlignedRead:
    def __init__(
        self,
        read_id,
        contig_id,
        sequence,
        read_start,
        read_stop,
        contig_start,
        contig_stop,
    ):
        self.read_id = read_id
        self.contig_id = contig_id
        self.sequence = sequence
        self.read_start = read_start
        self.read_stop = read_stop
        self.contig_start = contig_start
        self.contig_stop = contig_stop


def get_contig_start(kmer, kmers, contig_start):
    overlap = 0
    while kmers[kmer - overlap - 1].sufix == kmers[kmer].prefix:
        overlap += 1
    return contig_start + overlap


def assemble_contigs(contig_kmers, read_dict):
    contigs = []
    for contig_id, kmers in contig_kmers.items():
        aligned_reads = []
        sequence = kmers[0].sequence
        read_id = kmers[0].read_id
        read_start = kmers[0].start
        contig_start = 0
        for kmer in range(1, len(kmers)):
            read_num = kmers[kmer].read_id
            if read_num != read_id:
                read_stop = kmers[kmer - 1].stop
                # contig_start = get_contig_start(kmer-1, kmers, contig_start)
                contig_stop = contig_start + len(sequence)
                sequence = read_dict[read_id][read_start:read_stop]
                read = AlignedRead(
                    read_id=read_id,
                    contig_id=contig_id,
                    sequence=sequence,
                    read_start=read_start,
                    read_stop=read_stop,
                    contig_start=contig_start,
                    contig_stop=contig_stop,
                )
                read_id = read_num
                read_start = kmers[kmer].start
                contig_start = contig_stop + 1
                aligned_reads.append(read)
        sequence += kmers[kmer - 1].sequence[-1]
        contig = Contig(
            contig_id=contig_id,
            aligned_reads=aligned_reads,
            sequence=sequence,
            direction=1,
        )
        contigs.append(contig)
    return contigs


"""
def get_reverse_contigs(contigs, read_dict):
    rvs_contigs = []
    contig_id = len(contigs)
    for contig in contigs:
        sequence = contig.sequence[::-1]
        aligned_reads = []
        for read in contig.aligned_reads:
            read_id = read.read_id
            read_start = read.read_stop
            read_stop = read.read_start
            sequence = read_dict[read_id][read_start:read_stop]
            contig_start = read.contig_start
            contig_stop = read.contig_stop
            read = AlignedRead(
                    read_id = read_id, 
                    contig_id=contig_id, 
                    sequence = sequence,
                    read_start = read_start, 
                    read_stop = read_stop, 
                    contig_start = contig_start, 
                    contig_stop = contig_stop
                    )
        read_id = read_num
        rvs_contig = Contig(contig_id = contig_id, aligned_reads = aligned_reads, sequence = sequence, direction = -1)
        rvs_contigs.append(rvs_contig)
        contig_id += 1
    return(rvs_contigs) 

"""


def assembly(path_dict, read_kmers, read_dict):
    contig_kmers = get_contig_kmers(path_dict, read_kmers)
    contigs = assemble_contigs(contig_kmers, read_dict)
    # rvs_contigs = get_reverse_contigs(contigs, read_dict)
    return contigs
