__author__ = "Keenan Manpearl"
__date__ = "2023/04/09"

"""
Functions needed to format output files    
"""

import pandas as pd


def get_longest_contig(contigs, out_dir):
    length_longest_contig = 0
    for contig in contigs:
        if len(contig.sequence) > length_longest_contig:
            length_longest_contig = len(contig.sequence)
            longest_contig_id = contig.contig_id
            longest_contig_sequence = contig.sequence
            longest_contig = contig
    with open(f"{out_dir}/ALLELES.fasta", "w") as file:
        file.write(f">contig_{longest_contig_id}\n{longest_contig_sequence}")
    return longest_contig


def format_aln_file(longest_contig, out_dir):
    reads = longest_contig.aligned_reads
    aln = pd.DataFrame(reads).T
    aln.to_csv(f"{out_dir}/ALLELES.aln", index=False)


def save_required_ouputs(contigs, out_dir):
    longest_contig = get_longest_contig(contigs, out_dir)
    format_aln_file(longest_contig, out_dir)
