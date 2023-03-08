__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

import argparse as arg 
import pandas as pd

from data_loader import parse_query, parse_reads
from alignment import get_reads_to_align, get_reads_to_assemble
from assembly import create_graph, create_adjacency_matrix
from kmers import get_kmers


parser = arg.ArgumentParser(description='Align sequence reads to a query sequence and assemble contigs from aligned reads')
required = parser.add_argument_group('required arguments')
required.add_argument('--q', '-query_file', type=str, help='path to the query FASTA file', required = True)
required.add_argument('--r', '-read_file', type=str, help='path to the reads FASTA file', required = True)
parser.add_argument('--k', '-kmer_length', type=str, help='length of kmers', default = 10)
parser.add_argument('--m', '-match_score', type=bool, help='alignment score for matching base pairs', default = 1)
parser.add_argument('--mi', '-mismatch_score', type=bool, help='alignment score for non-matching base pairs', default = -1)
parser.add_argument('--g', '-gap_score', type=bool, help='alignment score for introducing a gap', default = -1)
parser.add_argument('--t', '-score_threshold', type=bool, help='minimum normalized score needed to be considered an alignment', default = 0.75)
parser.add_argument('--s', '-save', type=bool, help='if True, save a csv file of the graph', default = False)
args = parser.parse_args()


query_seq = parse_query(args.q)
read_dict = parse_reads(args.r)
k = args.k
match_score = args.m 
gap_score = args.g 
mismatch_score = args.mi 
threshold = args.t 
save = args.s


'''
# for testing purposes 
query_seq = parse_query("test_data/query.FASTA.txt")
read_dict = parse_reads("test_data/reads.FASTA.txt")
k = 10
match_score = 1 
gap_score = -1
mismatch_score = -1 
threshold = .75 
save = True
'''

# TODO: new function, modify to include reverse reads 
reads_to_align, no_alignment = get_reads_to_align(query_seq, read_dict)
fwd_reads_to_assemble, reverse_reads_to_assemble, no_alignment = get_reads_to_assemble(query_seq, reads_to_align, no_alignment, match_score, gap_score, mismatch_score, threshold, save)
fwd_kmers = get_kmers(fwd_reads_to_assemble, k)
reverse_kmers = get_kmers(reverse_reads_to_assemble, k)
# TODO: add reverse kmers to graph
graph = create_graph(fwd_kmers)
if save == True:
    graph_df = pd.DataFrame(graph, columns = ["source", "target"]) 
    graph_df.to_csv('./output/test_graph.csv')

