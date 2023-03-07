__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

import numpy as np
from kmers import create_query_kmers, create_reads_kmers


def get_reads_to_align(query_seq, read_dict,k=3):
    query_kmers = create_query_kmers(k, query_seq)
    read_kmers = create_reads_kmers(k,read_dict)
    reads = []
    for read_k in read_kmers:
        for query_k in query_kmers:
            if read_k == query_k:
                reads.append(read_k.read_id)
    reads_to_align = {key: read_dict[key] for key in np.unique(reads)}
    return(reads_to_align)



def score(base1, base2, match_score, mismatch_score):
    if base1 == base2:
        return match_score
    else:
        return mismatch_score

def compare_sequences(query, sequence, match_score=10, gap_score=-5,mismatch_score=-5, threshold=50):
    no_alignment = False
    scores = np.zeros((len(query)+1, len(sequence)+1))
    best_score = 0
    best_score_idx = (0,0)
    for i in range(1, len(query)+1):
        for j in range(1, len(sequence)+1):
            scores[i][j] = max(
                scores[i][j-1] + gap_score,
                scores[i-1][j] + gap_score,
                scores[i-1][j-1] + score(query[i-1], sequence[j-1], match_score, mismatch_score),
                0
            )
            if scores[i][j] >= best_score: 
                best_score = scores[i][j] 
                best_score_idx = [i,j]
    if best_score < threshold:
        scores = None
        no_alignment = True
    elif best_score_idx[0] != len(query) and best_score_idx[1] != len(sequence):
        scores = None
        no_alignment = True
    return scores, no_alignment


def get_reads_to_assemble(query, read_dict):
    score_list = {}
    score_list_reverse = {}
    no_alignment = []
    fwd_reads_to_assemble = {}
    reverse_reads_to_assemble = {}
    for read, sequence in read_dict.items():
        scores, no_alignment = compare_sequences(query, sequence)
        score_list[read] = scores 
        reverse = sequence[::-1]
        scores, no_alignment_reverse = compare_sequences(query, reverse)
        score_list_reverse[read] = scores
        if no_alignment:
            fwd_reads_to_assemble[read] = sequence 
        if no_alignment_reverse:
            reverse_reads_to_assemble[read] = sequence
    if len(fwd_reads_to_assemble) == 0 and len(reverse_reads_to_assemble) == 0:
        raise Exception("No reads align to query sequence")
    else:
        return fwd_reads_to_assemble, reverse_reads_to_assemble
