__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

"""
Classes and functions needed to create kmers
"""


class Kmer:
    """
    classes to hold information about kmers
    """

    def __init__(self, read_id, id, sequence, prefix, suffix, start, stop, direction):
        self.read_id = read_id
        self.id = id
        self.sequence = sequence
        self.prefix = prefix
        self.sufix = suffix
        self.start = start
        self.stop = stop
        self.direction = direction


def create_reads_kmers(read_dict, k, id, direction):
    """
    create kmers for reads

    Args:
        read_dict (dict): read ids (keys) and sequences (values)
        k (int): length of kmers
        id (int): number to use as starting id for first kmer (0 for forward direction)
        direction (int): 1 for forward or -1 for reverse

    Returns:
        list : instances of class Kmer made from reads, and id to be used for start id for reverse kmers
    """
    kmers = []
    for read_id, read in read_dict.items():
        for i in range(len(read) - k + 1):
            start = i
            stop = i + k
            sequence = read[start:stop]
            prefix = sequence[:-1]
            sufix = sequence[1:]
            kmer = Kmer(read_id, id, sequence, prefix, sufix, start, stop, direction)
            kmers.append(kmer)
            id += 1
    return kmers, id


def create_query_kmers(query, k):
    """
    create kmers from query sequence

    Args:
        query (str): sequence of query
        k (int): size of kmer

    Returns:
        list: instances of class Kmer that make up query sequence
    """
    kmers = []
    id = 0
    for i in range(len(query) - k + 1):
        start = i
        stop = i + k
        sequence = query[start:stop]
        prefix = sequence[:-1]
        sufix = sequence[1:]
        kmer = Kmer("query", id, sequence, prefix, sufix, start, stop, 1)
        kmers.append(kmer)
        id += 1
    return kmers


def get_all_kmers(read_dict, rvs_read_dict, k):
    """
    generate forward and reverse read kmers

    Args:
        read_dict (dict): forward reads with read ids (key) and sequences (value)
        rvs_read_dict (dict): reverse reads with read ids (key) and sequences (value)
        k (int): size of kmers

    Returns:
        list : instances of class Kmers that make up all forward and reverse reads
    """
    fwd_kmers, id = create_reads_kmers(read_dict, k, 0, 1)
    # rvs_kmers, id = create_reads_kmers(rvs_read_dict, k, id, -1)
    return fwd_kmers  # + rvs_kmers


def create_contig_kmers(contigs, k):
    """
    create kmers of contigs

    Args:
        contigs (list): instances of class Contigs assembled from reads
        k (int): size of kmers

    Returns:
        list: instances of class Kmer assembled from all contigs
    """
    kmers = []
    id = 0
    for contig in contigs:
        contig_sequence = contig.sequence
        for i in range(len(contig_sequence) - k + 1):
            start = i
            stop = i + k
            sequence = contig_sequence[start:stop]
            prefix = sequence[:-1]
            sufix = sequence[1:]
            kmer = Kmer(contig.contig_id, id, sequence, prefix, sufix, start, stop, 1)
            kmers.append(kmer)
            id += 1
    return kmers
