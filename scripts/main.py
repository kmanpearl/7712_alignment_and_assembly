__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"


import pandas as pd

from data_loader import parse_query, parse_reads
from alignment import get_reads_to_assemble
from assembly import create_graph
from kmers import Kmer, create_reads_kmers

query_seq = parse_query("./test_data/query.FASTA.txt")
read_dict = parse_reads("./test_data/reads.FASTA.txt")



# TODO: figure out why this function currently misses some reads that should align
# reads_to_align = get_reads_to_align(query_seq, read_dict)
fwd_reads_to_assemble, reverse_reads_to_assemble = get_reads_to_assemble(query_seq, read_dict)
kmers = create_reads_kmers(4, fwd_reads_to_assemble)
for kmer in kmers:
    Kmer.define_prefix_and_sufix(kmer)

graph = create_graph(kmers)
graph_df = pd.DataFrame(graph, columns = ["source", "target"])
# uncomment to save graph for visualization purposes 
graph_df.to_csv('./output/graph.csv')

