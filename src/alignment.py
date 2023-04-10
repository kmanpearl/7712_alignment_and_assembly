__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

"""
Functions needed to align reads to query sequence 
"""

import csv
import datetime

import numpy as np

from src.kmers import create_contig_kmers, create_query_kmers


def get_contigs_to_align(query_seq, contigs, k):
    """
    find all contigs that have a kmer match to the query sequence

    Args:
        query_seq (str): sequence to align against
        contigs (list): instances of class Contigs containing all assemble reads
        k (int): kmer size

    Raises:
        Exception: no kmer matches

    Returns:
        touple : list of instances of class Contig that contain a kmer match, list of class Contigs that do not share a kmer match
    """
    query_kmers = create_query_kmers(query_seq, k)
    contig_kmers = create_contig_kmers(contigs, k)
    kmer_match = []
    no_alignment = []
    contigs_to_align = []
    for contig_kmer in contig_kmers:
        for query_kmer in query_kmers:
            if contig_kmer.sequence == query_kmer.sequence:
                contig_id = contig_kmer.read_id
                if contig_id not in kmer_match:
                    kmer_match.append(contig_id)
                if contig_id in no_alignment:
                    no_alignment.remove(contig_id)
                break
    for contig in contigs:
        if contig.contig_id in kmer_match:
            contigs_to_align.append(contig)
    if len(contigs_to_align) == 0:
        raise Exception("No contigs align to query sequence")
    print(
        f"{datetime.datetime.now()}: found {len(contigs_to_align)} contigs for alignment"
    )
    return contigs_to_align, no_alignment


def score_matches(base1, base2, match_score, mismatch_score):
    """
    determine if two base pairs match and return the proper score

    Args:
        base1 (str): base pair
        base2 (str): base pair
        match_score (_type_): _description_
        mismatch_score (_type_): _description_

    Returns:
        int: score for the alignment
    """
    if base1 == base2:
        return match_score
    else:
        return mismatch_score


def compare_sequences(
    query, sequence, match_score, gap_score, mismatch_score, threshold
):
    """
    generate an alignmnet score between two sequences

    Args:
        query (str): query sequence to align against
        sequence (str): contig to align
        match_score (int): score for matching bases
        gap_score (int): penalty for introducing a gap
        mismatch_score (int): penalty for a mismatch
        threshold (float): minimum score to be considered an alignmnet (between 0-1)

    Returns:
        touple: float of alignment score, bool of whether alignmnet was above the threshold
    """
    alignment = True
    scores = np.zeros((len(query) + 1, len(sequence) + 1))
    best_score = 0
    # best_score_idx = (0, 0)
    for i in range(1, len(query) + 1):
        for j in range(1, len(sequence) + 1):
            max_score = max(
                0,
                scores[i - 1][j - 1]
                + score_matches(
                    query[i - 1], sequence[j - 1], match_score, mismatch_score
                ),
                scores[i - 1][j] + gap_score,
                scores[i][j - 1] + gap_score,
            )
            scores[i][j] = max_score
            if scores[i][j] >= best_score:
                best_score = scores[i][j]
                # best_score_idx = (i, j)
    # if best read score normalized for sequence length is below threshold
    # add to list of unaligned reads.
    # or return score
    score = best_score / len(sequence)
    if score < threshold:
        best_score = None
        alignment = False
        return best_score, alignment
    else:
        return score, alignment


def alignment(
    query_seq, contigs, match_score, gap_score, mismatch_score, threshold, save
):
    """
    align all contigs against the query

    Args:
        query_seq (str): query sequence to align against
        contigs (list): instances of class Contig containing all assembled contigs
        match_score (int): score for matching bases
        gap_score (int): penalty for introducing a gap
        mismatch_score (int): penalty for a mismatch
        threshold (float): minimum score to be considered an alignmnet (between 0-1)
        save (bool): whether to save csv of all alignment scores (for true alignments)

    Raises:
        Exception: _description_

    Returns:
        _type_: _description_
    """
    no_alignment = []
    score_dict = {}
    aligned_contigs = []
    contigs_to_align, no_alignment = get_contigs_to_align(query_seq, contigs, k=5)
    num_contigs = len(contigs_to_align)
    for i in range(num_contigs):
        contig_id = contigs_to_align[i].contig_id
        sequence = contigs_to_align[i].sequence
        score, alignment = compare_sequences(
            query_seq, sequence, match_score, gap_score, mismatch_score, threshold
        )
        if alignment:
            score_dict[contig_id] = score
            aligned_contigs.append(contigs_to_align[i])
            if contig_id in no_alignment:
                no_alignment.remove(contig_id)
    if len(aligned_contigs) == 0:
        raise Exception("No contigs align to query sequence")
    if save:
        with open("output/alignment_scores.csv", "w") as csv_file:
            writer = csv.writer(csv_file)
            for key, value in score_dict.items():
                writer.writerow([key, value])
    return aligned_contigs, no_alignment
