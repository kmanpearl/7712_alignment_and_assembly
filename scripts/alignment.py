__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

import csv
import numpy as np
import pandas as pd
from kmers import create_query_kmers, create_reads_kmers


def get_reads_to_align(query_seq, read_dict, k=5):
    query_kmers = create_query_kmers(k, query_seq)
    read_kmers = create_reads_kmers(k,read_dict)
    aligned = []
    no_alignment = []
    for read_k in read_kmers:
        for query_k in query_kmers:
            if read_k == query_k:
                aligned.append(read_k.read_id)
                break
    aligned = np.unique(aligned)
    reads_to_align = {key: read_dict[key] for key in aligned}
    all_reads = list(read_dict.keys())
    for read in all_reads:
        if read not in aligned:
            no_alignment.append(read)
    return reads_to_align, no_alignment



def score(base1, base2, match_score, mismatch_score):
    if base1 == base2:
        return match_score
    else:
        return mismatch_score

def compare_sequences(query, sequence, match_score, gap_score, mismatch_score, threshold):
    alignment = True
    scores = np.zeros((len(query)+1, len(sequence)+1))
    best_score = 0
    best_score_idx = (0,0)
    for i in range(1, len(query)+1):
        for j in range(1, len(sequence)+1):
            max_score = max(
                0,
                scores[i-1][j-1] + score(query[i-1], sequence[j-1], match_score, mismatch_score),
                scores[i-1][j] + gap_score,
                scores[i][j-1] + gap_score,
            )
            scores[i][j] = max_score
            if scores[i][j] >= best_score: 
                best_score = scores[i][j] 
                best_score_idx = (i,j)
    # if best read score normalized for sequence length is below threshold
    # add to list of unaligned reads. 
    if best_score/len(sequence) < threshold:
        scores = None
        alignment = False
    # the end of the read must align to somewhere in the query 
    # or anywhere in the query must align to the end of the read 
    # if only part of the read aligns to the middle of the query it is not a true alignment 
    elif best_score_idx[0] != len(query) and best_score_idx[1] != len(sequence):
        scores = None
        alignment = False
    return best_score/len(sequence), alignment


def get_reads_to_assemble(query, read_dict, no_alignment, match_score, gap_score, mismatch_score, threshold, save):
    fwd_score_dict = {}
    reverse_score_dict = {}
    fwd_reads_to_assemble = {}
    reverse_reads_to_assemble = {}
    for read, sequence in read_dict.items():
        score, alignment = compare_sequences(query, sequence, match_score, gap_score, mismatch_score, threshold)
        fwd_score_dict[read] = score
        reverse = sequence[::-1]
        reverse_score, alignment_reverse = compare_sequences(query, reverse, match_score, gap_score, mismatch_score, threshold)
        reverse_score_dict[read] = reverse_score
        if alignment and score > reverse_score:
            fwd_reads_to_assemble[read] = sequence 
        if alignment_reverse and reverse_score > score:
            reverse_reads_to_assemble[read] = sequence[::1]
        if not alignment and not alignment_reverse:
            no_alignment.append(read)
    if len(fwd_reads_to_assemble) == 0 and len(reverse_reads_to_assemble) == 0:
        raise Exception("No reads align to query sequence")

    if save == True:
        # TODO: updated csv file to include absolute and normalized scores
        with open('./output/test_fwd_alignment_scores.csv', 'w') as csv_file:  
            writer = csv.writer(csv_file)
            for key, value in fwd_score_dict.items():
                writer.writerow([key, value])
        with open('./output/test_reverse_alignment_scores.csv', 'w') as csv_file:  
            writer = csv.writer(csv_file)
            for key, value in reverse_score_dict.items():
                writer.writerow([key, value])
    return fwd_reads_to_assemble, reverse_reads_to_assemble, no_alignment
