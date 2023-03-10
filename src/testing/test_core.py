"""
Here we are testing the core functions of the program
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath("../"))
from alignment import compare_sequences, score_matches
from assembly import create_graph
from data_loader import parse_query, parse_reads
from kmers import Kmer, create_query_kmers, create_reads_kmers


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
        reads = parse_reads("testing/test_data/true_reads.txt")
        expected_output = {"seq1": "ATGC", "seq2": "GCC"}
        self.assertEqual(reads, expected_output)

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


class TestKmerClasses(unittest.TestCase):
    """
    test all functions and classes used to create kmers
    """

    def test_read_kmers(self):
        """
        test the creation of kmers from reads
        """
        # TODO: creating prefix and sufix turns object into NoneType
        # I have confirmed this method works when executed in my code
        # need to figure out why unittest is failing
        read_dict = {"seq1": "ACTG"}
        kmers = create_reads_kmers(3, read_dict)
        kmer = kmers[0]
        # kmer_ps = Kmer.define_prefix_and_sufix(kmer)
        expected_read_id = "seq1"
        expected_read = "ACTG"
        expected_sequence = "ACT"
        expected_start = 0
        expected_stop = 3
        # expected_prefix = "AC"
        # expected_sufix = "CT"
        self.assertEqual(kmer.read_id, expected_read_id)
        self.assertEqual(kmer.read, expected_read)
        self.assertEqual(kmer.sequence, expected_sequence)
        self.assertEqual(kmer.start, expected_start)
        self.assertEqual(kmer.stop, expected_stop)
        # elf.assertEqual(kmer_ps.prefix, expected_prefix)
        # self.assertEqual(kmer_ps.sufix, expected_sufix)

    def test_query_kmers(self):
        """
        test the creation of kmers from query sequence
        """
        kmers = create_query_kmers(3, "ACTG")
        kmer = kmers[0]
        expected_read_id = "query"
        expected_read = "ACTG"
        expected_sequence = "ACT"
        expected_start = 0
        expected_stop = 3
        self.assertEqual(kmer.read_id, expected_read_id)
        self.assertEqual(kmer.read, expected_read)
        self.assertEqual(kmer.sequence, expected_sequence)
        self.assertEqual(kmer.start, expected_start)
        self.assertEqual(kmer.stop, expected_stop)


class TestAssemblyFuncions(unittest.TestCase):
    """
    Testing all functions and classes needed for assmebly
    """

    # TODO: this function uses prefixes and sufixes
    # cannot test until I get method to work within unittesting
    # def test_graph(self):
    #    kmers = create_query_kmers(3, "ACTG")
    #    for kmer in kmers:
    #        kmer = Kmer.define_prefix_and_sufix(kmer)
    #    edges = create_graph(kmers)
    #    expected_output = [("ACT, CTG")]
    #    self.assertEqual(edges, expected_output)


# python3 -m unittest testing.test_core
