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
import datetime

from src.alignment import alignment
from src.assembly import assembly, graph_traversal
from src.data_loader import parse_query, parse_reads
from src.kmers import get_all_kmers
from src.output import save_required_ouputs

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
    "--k",
    "-kmer_size",
    type=int,
    help="length of kmers: must be shorter than the shortest read",
    default=15,
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
    "--t",
    "-score_threshold",
    type=float,
    help="minimum normalized score needed to be considered an alignment: value must be between 0-1",
    default=0.5,
)
parser.add_argument(
    "--s", "-save", type=bool, help="if True, save intermediate outputs", default=False
)
args = parser.parse_args()


query_seq = parse_query(args.q)
print(f"{datetime.datetime.now()}: parsed query files")
read_dict, rvs_read_dict = parse_reads(args.r)
print(f"{datetime.datetime.now()}: parsed reads files")
k = args.k
match_score = args.m
gap_score = args.g
mismatch_score = args.mi
threshold = args.t
save = args.s


read_kmers = get_all_kmers(read_dict, rvs_read_dict, k)
print(f"{datetime.datetime.now()}: created kmers")
path_dict = graph_traversal(read_kmers, save)
all_contigs = assembly(path_dict, read_kmers, read_dict)
aligned_contigs, no_alignment = alignment(
    query_seq=query_seq,
    contigs=all_contigs,
    match_score=match_score,
    gap_score=gap_score,
    mismatch_score=mismatch_score,
    threshold=threshold,
    save=save,
)
print(f"{datetime.datetime.now()}: formatting output")
save_required_ouputs(aligned_contigs)
print(f"{datetime.datetime.now()}: done")
