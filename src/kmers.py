__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

"""
Classes and functions needed to create kmers
"""


class Kmer:
    def __init__(self, read_id, id, sequence, prefix, suffix, start, stop, direction):
        self.read_id = read_id
        self.id = id
        self.sequence = sequence
        self.prefix = prefix
        self.sufix = suffix
        self.start = start
        self.stop = stop
        self.direction = direction

    # def __eq__(self, other):
    #    return self.sequence == other.sequence


def create_reads_kmers(read_dict, k, id, direction):
    """
     Parameters
    ----------
    k : INT
        size of kmer
    read_dict : DICT
        read_id (keys) and sequence (values)

    Returns
    -------
    kmers : LIST
        list of Kmer objects
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


"""
def get_kmers(read_dict, k):
    kmers = create_reads_kmers(k, read_dict)
    return kmers
"""


def get_all_kmers(read_dict, rvs_read_dict, k):
    fwd_kmers, id = create_reads_kmers(read_dict, k, 0, 1)
    rvs_kmers, id = create_reads_kmers(rvs_read_dict, k, id, -1)
    return fwd_kmers + rvs_kmers


def create_contig_kmers(contigs, k):
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
