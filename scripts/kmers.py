__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

class Kmer:
    def __init__(self, read_id, read, sequence, start, stop):
        self.read_id = read_id
        self.read = read
        self.sequence = sequence
        # TODO: do i need these?
        self.start = start
        self.stop = stop  
    
    def __eq__(self, other):
        return self.sequence == other.sequence

    def define_prefix_and_sufix(self):
        self.prefix = self.sequence[:-1]
        self.sufix = self.sequence[1:] 

def create_reads_kmers(k, read_dict):
    '''
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
    '''
    kmers = []
    for read_id, read in read_dict.items():
        for i in range(len(read)-k+1):
            start = i
            stop = i+k
            sequence = read[start:stop]
            kmer = Kmer(read_id, read, sequence, start, stop)
            kmers.append(kmer)
    return(kmers)

def create_query_kmers(k, query):
    kmers = []
    for i in range(len(query)-k+1):
        start = i
        stop = i+k
        sequence = query[start:stop]
        kmer = Kmer("query", query, sequence, start, stop)
        kmers.append(kmer)
    return(kmers)

def get_kmers(reads_to_assemble, k):
    if len(reads_to_assemble) > 0:
        kmers = create_reads_kmers(k, reads_to_assemble)
        for kmer in kmers:
            Kmer.define_prefix_and_sufix(kmer)
        return(kmers)
    else:
        return None 
