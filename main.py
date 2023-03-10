__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

"""
This script takes a query sequence and collection of DNA sequence reads.
First sequence reads and the query sequence are split into kmers. 
If any sequence has an exact kmer match to a query kmer, 
it is aligned agains the query using dynamic programming.
Then any reads with an alignmnet score above a user inputed threshold 
are assembled into contigs using a de Bruijn graph. 
"""

import argparse as arg

import pandas as pd

from src.alignment import get_reads_to_align, get_reads_to_assemble
from src.assembly import create_graph
from src.data_loader import parse_query, parse_reads
from src.kmers import get_kmers

parser = arg.ArgumentParser(
    description="Align sequence reads to a query sequence and assemble contigs from aligned reads"
)
required = parser.add_argument_group("required arguments")
required.add_argument(
    "--q", "-query_file", type=str, help="path to the query FASTA file", required=True
)
required.add_argument(
    "--r", "-read_file", type=str, help="path to the reads FASTA file", required=True
)
parser.add_argument(
    "--ka",
    "-alignment_kmer",
    type=int,
    help="length of kmers to use for alignment",
    default=5,
)
parser.add_argument(
    "--kb",
    "-assembly_kmer",
    type=int,
    help="length of kmers to use for assembly",
    default=10,
)
parser.add_argument(
    "--m",
    "-match_score",
    type=int,
    help="alignment score for matching base pairs",
    default=1,
)
parser.add_argument(
    "--mi",
    "-mismatch_score",
    type=int,
    help="alignment score for non-matching base pairs",
    default=-1,
)
parser.add_argument(
    "--g",
    "-gap_score",
    type=int,
    help="alignment score for introducing a gap",
    default=-1,
)
parser.add_argument(
    "--st",
    "-score_threshold",
    type=float,
    help="minimum normalized score needed to be considered an alignment",
    default=0.75,
)
parser.add_argument(
    "--s", "-save", type=bool, help="if True, save intermediate outputs", default=False
)
args = parser.parse_args()


query_seq = parse_query(args.q)
read_dict = parse_reads(args.r)
k_alignment = args.ka
k_assembly = args.kb
match_score = args.m
gap_score = args.g
mismatch_score = args.mi
threshold = args.st
save = args.s


"""
# for testing purposes
query_seq = parse_query("sample_data/query.FASTA.txt")
read_dict = parse_reads("sample_data/reads.FASTA.txt")
k_alignment = 3
k_assembly = 5
match_score = 1
gap_score = -1
mismatch_score = -1
threshold = 0.75
save = True
"""
# TODO: new function, modify to include reverse reads
reads_to_align, no_alignment = get_reads_to_align(query_seq, read_dict, k_assembly)
# no_alignment = []
fwd_reads_to_assemble, reverse_reads_to_assemble, no_alignment = get_reads_to_assemble(
    query_seq,
    reads_to_align,
    no_alignment,
    match_score,
    gap_score,
    mismatch_score,
    threshold,
    save,
)
fwd_kmers = get_kmers(fwd_reads_to_assemble, k_alignment)
reverse_kmers = get_kmers(reverse_reads_to_assemble, k_alignment)
# TODO: add reverse kmers to graph
graph = create_graph(fwd_kmers)
if save == True:
    graph_df = pd.DataFrame(graph, columns=["source", "target"])
    graph_df.to_csv("output/graph.csv")
