__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

"""
Functions needed to read in input files    
"""


def parse_query(fp):
    """
    parse query sequence from fasta file

    Args:
        fp (str): path to query fasta

    Raises:
        Exception: is less than two lines
        Exception: if more than two lines
        Exception: if file contains any letters besided ACTG
        Exception: if sequence id lines do not start with '>'

    Returns:
        str: sequence of query
    """
    nucleotides = set("ATCG")
    with open(fp, "r") as file:
        lines = file.readlines()
        if len(lines) < 2:
            raise Exception("Query sequence must be a FASTA file")
        elif len(lines) > 2:
            raise Exception("Query sequence must be a FASTA file with only one entry")
        if lines[0].startswith(">"):
            sequence = lines[1].upper().strip()
            if set(sequence) > nucleotides:
                raise Exception("Query must be a DNA sequence")
            return sequence
        else:
            raise Exception("Query sequence must be a FASTA file")


def parse_reads(fp):
    """
    create dictionaries of reads

    Args:
        fp (str): path to reads fasta file

    Raises:
        Exception: if unequal number of lines
        Exception: if only one read
        Exception: if sequence id lines do not start with '>'
        Exception: if file contains any letters besided ACTG

    Returns:
        dicts: dictionaries of forward and reverse reads with read id (key) and sequence (value)
    """
    nucleotides = set("ATCG")
    fwd_read_dict = {}
    rvs_read_dict = {}
    with open(fp, "r") as file:
        lines = file.readlines()
        if len(lines) % 2 != 0:
            raise Exception("Reads must be a FASTA file")
        elif len(lines) == 2:
            raise Exception(
                "Reads FASTA file must contain more than one sequence to assemble"
            )
        i = 0
        while i < len(lines):
            read_id = lines[i]
            if not read_id.startswith(">"):
                raise Exception("Reads must be a FASTA file")
            # remove > from line and any trailing whitespace
            read_id = read_id[1:].strip()
            sequence = lines[i + 1].upper().strip()
            if set(sequence) > nucleotides:
                raise Exception("Reads must be DNA sequences")
            fwd_read_dict[read_id] = sequence
            rvs_read_dict[read_id] = sequence[::-1]
            i += 2

    return fwd_read_dict, rvs_read_dict
