__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

"""
Functions needed to assemble sequence reads 
"""


class Contig:
    """
    holds information about assembled contigs
    """

    def __init__(self, contig_id, aligned_reads, sequence, direction):
        self.contig_id = contig_id
        self.sequence = sequence
        self.direction = direction
        self.aligned_reads = aligned_reads

    def add_kmers(self, kmers):
        self.kmers = kmers


def get_contig_kmers(all_paths, read_kmers):
    """
    get the id of all kmers that make up contigs

    Args:
        all_paths (list): instances of class Path that represent all possible paths through the graph
        read_kmers (list): instances of class Kmer that make up all reads

    Returns:
        dict: contig ids (keys) and ids of all kmers that make up contig (values)"""
    contig_kmers = {}
    for path in all_paths:
        contig_id = path.contig_id
        node_kmers = {}
        for node in path.path:
            kmers = []
            for kmer in read_kmers:
                if kmer.sequence == node:
                    kmers.append(kmer)
            kmers.sort(key=lambda x: x.id)
            node_kmers[node] = kmers
        if contig_id in contig_kmers.keys():
            contig_kmers[contig_id].append(node_kmers)
        else:
            contig_kmers[contig_id] = node_kmers
    return contig_kmers


def assemble_contigs(contig_kmers):
    contigs = []
    for contig_id, node_dict in contig_kmers.items():
        reads = {}
        sequence = str()
        contig_position = 0
        for node, kmers in node_dict.items():
            for kmer in kmers:
                read_id = kmer.read_id
                if read_id not in reads.keys():
                    reads[read_id] = {
                        "sseqid": read_id,
                        "qseqid": contig_id,
                        "sstart": kmer.start,
                        "send": kmer.stop,
                        "qstart": contig_position,
                        "qend": contig_position,
                    }
                else:
                    reads[read_id]["send"] = kmer.stop
                    reads[read_id]["qend"] = contig_position
                if len(sequence) == 0:
                    sequence = node[:-1]
                    contig_position = len(sequence)
            sequence += node[-1]
            contig_position += 1
        contig = Contig(
            contig_id=contig_id,
            aligned_reads=reads,
            sequence=sequence,
            direction=1,
        )
        contigs.append(contig)
    return contigs


def assembly(paths, read_kmers):
    """
    assemble all contigs

    Args:
        paths (list): instances of class Path that represent all possible paths through the graph
        read_kmers (list): instances of class Kmer that make up all reads

    Returns:
        list: instances of class Contigs
    """
    contig_kmers = get_contig_kmers(paths, read_kmers)
    contigs = assemble_contigs(contig_kmers)
    return contigs
