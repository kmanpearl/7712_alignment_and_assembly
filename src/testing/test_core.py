"""
Here we are testing the core functions of the program
"""

import os
import sys
import unittest

import pandas as pd

sys.path.append(os.path.abspath("../"))
from alignment import alignment, compare_sequences, get_contigs_to_align, score_matches
from assembly import Contig, get_contig_kmers
from data_loader import parse_query, parse_reads
from graph import (
    create_adjacency_matrix,
    create_graph,
    find_all_paths,
    find_start_stop_nodes,
)
from kmers import create_query_kmers, create_reads_kmers


class TestReaderFunctions(unittest.TestCase):
    """
    Testing related to the data_loader functions.
    Test that begin with test_reads are testing the parse_reads function
    Tests that begin with test_query are testing the parse_query function
    """

    def test_reads_true(self):
        """
        test that the expected ouput is returned
        """
        fwd_read_dict, rvs_read_dict = parse_reads("testing/test_data/true_reads.txt")
        expected_fwd_output = {"seq1": "ATGC", "seq2": "GCC"}
        expected_rvs_output = {"seq1": "CGTA", "seq2": "CCG"}
        self.assertEqual(fwd_read_dict, expected_fwd_output)
        self.assertEqual(rvs_read_dict, expected_rvs_output)

    def test_reads_length_fail(self):
        """
        test that an exception is raised if the file is not in fatsa format
        because it contains an uneven number of lines
        """
        self.assertRaises(
            Exception, parse_reads, "testing/test_data/read_length_fail.txt"
        )

    def test_reads_num_seq_fail(self):
        """
        test that an exception is raised if the file contains only one read
        """
        self.assertRaises(Exception, parse_reads, "testing/test_data/true_query.txt")

    def test_reads_ID_fail(self):
        """
        test that an exception is raise if the sequence ID does not start with >
        """
        self.assertRaises(Exception, parse_reads, "testing/test_data/reads_ID_fail.txt")

    def test_reads_DNA_fail(self):
        """
        test that an exception is raised if the file contains any letters
        other than ATCG
        """
        self.assertRaises(Exception, parse_reads, "testing/test_data/read_DNA_fail.txt")

    def test_query_true(self):
        """
        test that the expected output is returned
        """
        query = parse_query("testing/test_data/true_query.txt")
        expected_output = "ACTG"
        self.assertEqual(query, expected_output)

    def test_query_length_fail(self):
        """
        test that an exception is raised if the file conains < 2 lines
        """
        self.assertRaises(
            Exception, parse_query, "testing/test_data/query_length_fail.txt"
        )

    def test_query_num_seq_fail(self):
        """
        test that an exception is raised if the file contains > 2 lines
        """
        self.assertRaises(Exception, parse_query, "testing/test_data/true_reads.txt")

    def test_query_num_reads_fail(self):
        """
        test that an exception is raised if the file contains more than one read
        """
        self.assertRaises(Exception, parse_query, "testing/test_data/true_reads.txt")

    def test_query_DNA_fail(self):
        """
        test that an exception is raised if the file contains any letters
        other than ATCG
        """
        self.assertRaises(
            Exception, parse_query, "testing/test_data/query_DNA_fail.txt"
        )

    def test_query_ID_fail(self):
        """
        test that an exception is raised if the sequnce ID doesnt start with >
        """
        self.assertRaises(Exception, parse_query, "testing/test_data/query_ID_fail.txt")


class TestKmerClasses(unittest.TestCase):
    """
    test all functions and classes used to create kmers
    """

    def test_read_kmers(self):
        """
        test the creation of kmers from reads
        """
        read_dict = {"seq1": "ACTGAC"}
        kmers = create_reads_kmers(read_dict, 3, 0, 1)
        kmer_order = []
        for kmer in kmers[0]:
            kmer_order.append(kmer.sequence)
        kmer = kmers[0][0]
        expected_read_id = "seq1"
        expected_id = 0
        expected_sequence = "ACT"
        expected_start = 0
        expected_stop = 3
        expected_prefix = "AC"
        expected_sufix = "CT"
        expected_direction = 1
        expected_order = ["ACT", "CTG", "TGA", "GAC"]
        self.assertEqual(kmer.read_id, expected_read_id)
        self.assertEqual(kmer.id, expected_id)
        self.assertEqual(kmer.sequence, expected_sequence)
        self.assertEqual(kmer.start, expected_start)
        self.assertEqual(kmer.stop, expected_stop)
        self.assertEqual(kmer.prefix, expected_prefix)
        self.assertEqual(kmer.sufix, expected_sufix)
        self.assertEqual(kmer.direction, expected_direction)
        self.assertEqual(kmer_order, expected_order)

    def test_query_kmers(self):
        """
        test the creation of kmers from query sequence
        """
        kmers = create_query_kmers("ACTG", 3)
        kmer = kmers[0]
        expected_read_id = "query"
        expected_id = 0
        expected_sequence = "ACT"
        expected_start = 0
        expected_stop = 3
        expected_prefix = "AC"
        expected_sufix = "CT"
        expected_direction = 1
        self.assertEqual(kmer.read_id, expected_read_id)
        self.assertEqual(kmer.id, expected_id)
        self.assertEqual(kmer.sequence, expected_sequence)
        self.assertEqual(kmer.start, expected_start)
        self.assertEqual(kmer.stop, expected_stop)
        self.assertEqual(kmer.prefix, expected_prefix)
        self.assertEqual(kmer.sufix, expected_sufix)
        self.assertEqual(kmer.direction, expected_direction)


class TestReaderFunctions(unittest.TestCase):
    """
    testing of functions related to construction and traversal of de brujin graph
    """

    def test_graph(self):
        """
        test the creation of edge lists from list of kmers
        """

        kmers = create_reads_kmers({"seq1": "ACTGAC"}, 3, 0, 1)
        edges = create_graph(kmers[0])
        edge_kmers = []
        for edge in edges:
            edge_kmers.append((edge[0].sequence, edge[1].sequence))
        expected_output = [("ACT", "CTG"), ("CTG", "TGA"), ("TGA", "GAC")]
        self.assertEqual(edge_kmers, expected_output)

    def test_adjacency_matrix(self):
        """
        test the creation of an adjacency matrix from an edge list
        """
        kmers = create_reads_kmers({"seq1": "ACTGAC"}, 3, 0, 1)
        edges = create_graph(kmers[0])
        adj_matrix = create_adjacency_matrix(kmers[0], edges, False)
        expected_output = pd.DataFrame(
            {
                "ACT": [0, 0, 0, 0],
                "CTG": [1, 0, 0, 0],
                "TGA": [0, 1, 0, 0],
                "GAC": [0, 0, 1, 0],
            }
        )
        expected_output.index = ["ACT", "CTG", "TGA", "GAC"]
        self.assertEqual(adj_matrix.equals(expected_output), True)

    def test_start_stop_nodes(self):
        """
        test to find all possible start and stop nodes
        """
        adj_matrix = pd.DataFrame(
            {
                "ACT": [0, 0, 0, 0],
                "CTG": [1, 0, 0, 0],
                "TGA": [0, 1, 0, 0],
                "GAC": [0, 0, 1, 0],
            }
        )
        adj_matrix.index = ["ACT", "CTG", "TGA", "GAC"]
        start_nodes, stop_nodes = find_start_stop_nodes(adj_matrix)
        expected_start_node = ["ACT"]
        expected_stop_node = ["GAC"]
        self.assertEqual(start_nodes, expected_start_node)
        self.assertEqual(stop_nodes, expected_stop_node)

    def test_find_paths(self):
        """
        test to find paths through graph
        """
        adj_matrix = pd.DataFrame(
            {
                "ACT": [0, 0, 0, 0],
                "CTG": [1, 0, 0, 0],
                "TGA": [0, 1, 0, 0],
                "GAC": [0, 0, 1, 0],
            }
        )
        adj_matrix.index = ["ACT", "CTG", "TGA", "GAC"]
        paths = find_all_paths(adj_matrix)
        path_sequence = []
        for path in paths:
            path_sequence.append(path.path)
        expected_output = [["ACT", "CTG", "TGA", "GAC"]]
        self.assertEqual(path_sequence, expected_output)


class TestAssemblyFuncions(unittest.TestCase):
    """
    Testing all functions and classes needed for assmebly
    """

    def test_get_contig_kmers(self):
        """
        test to find all kmers in a contig
        """
        adj_matrix = pd.DataFrame(
            {
                "ACT": [0, 0, 0, 0],
                "CTG": [1, 0, 0, 0],
                "TGA": [0, 1, 0, 0],
                "GAC": [0, 0, 1, 0],
            }
        )
        adj_matrix.index = ["ACT", "CTG", "TGA", "GAC"]
        paths = find_all_paths(adj_matrix)
        kmers = create_reads_kmers({"seq1": "ACTGAC"}, 3, 0, 1)
        contig_kmers = get_contig_kmers(paths, kmers[0])
        kmer_ids = []
        for path in contig_kmers.values():
            for kmer in path:
                kmer_ids.append(kmer.id)
        expected_output = [0, 1, 2, 3]
        self.assertEqual(kmer_ids, expected_output)


class TestAlignmentFunctions(unittest.TestCase):
    """
    test all functions related to sequence alignment
    """

    def test_score_function_true(self):
        """
        check that the expected output is returned if bp are an exact match
        """
        score = score_matches("A", "A", 1, -1)
        expected_output = 1
        self.assertEqual(score, expected_output)

    def test_score_function_false(self):
        """
        check that the expected output is returned if bp are not an exact match
        """
        score = score_matches("A", "C", 1, -1)
        expected_output = -1
        self.assertEqual(score, expected_output)

    def test_compare_sequences_true(self):
        """
        test that two exact matches return a score of 1
        """
        alignment = compare_sequences("ACTG", "ACTG", 1, -1, -1, 0)
        expected_output = (1, True)
        self.assertEqual(alignment, expected_output)

    def test_compare_sequences_false(self):
        """
        test that two exact matches return a score of 0
        """
        alignment = compare_sequences("ACTG", "CGAT", 1, -1, -1, 0.5)
        expected_output = (None, False)
        self.assertEqual(alignment, expected_output)

    def test_get_contigs_to_align(self):
        """
        test function to get sequences to align based on kmer matches to the query
        """
        contig = [Contig(1, [0, 1, 2], "ACTGAC", 1)]
        contigs_to_align, no_alignment = get_contigs_to_align("ACTGAC", contig, 4)
        expected_output = contig
        self.assertEqual(contigs_to_align, expected_output)

    def test_no_contigs_to_align(self):
        """
        test that an exception is raised if their are no reads with kmer matches to query
        """
        contig = [Contig(1, [0, 1, 2], "ACTGAC", 1)]
        self.assertRaises(Exception, parse_reads, "TACGA", contig, 4)

    def test_alignment(self):
        """
        test that matching reads align to the query
        """
        contig = [Contig(1, [0, 1, 2], "ACTGAC", 1)]
        aligned_reads, no_alignment = alignment("ACTGAC", contig, 1, -1, -1, 0.5, False)
        expected_output = contig
        self.assertEqual(aligned_reads, expected_output)

    def test_alignment_false(self):
        """
        test that an exception is raised if no reads align to query
        """
        contig = [Contig(1, [0, 1, 2], "TACGA", 1)]
        self.assertRaises(Exception, alignment, "ACTGAC", contig, 1, -1, -1, 0.5, False)


# python3 -m unittest testing.test_core
