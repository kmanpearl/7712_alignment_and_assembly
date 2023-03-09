"""
Here we are testing the core functions of the program
"""

import unittest  # this imports the unit testing stuff

# source import
from src.alignment import get_reads_to_assemble


# next we need to create a class that will house all the tests
# -- we need to inheret the TestCase class in order to do our custom test
class TestCore(unittest.TestCase):
    # each method should have a `test_` prefix
    def test_reader(self):
        """
        - description of tests (what is testing)
        - what is the expected function
        - is this a postive, or a negative test()
            - postive = the expected is captured
                - Expected output is captured
                    - Example: the expected kmer matches with the generated kmer
            - negative = the expected incorrect is captured
                - expcetion handling (making things fail on purpose)
        """
        pass

    def test_2(self):
        """
        example of a postive test
        """
        # part1: your implementation output (test case)
        test_my_function = get_reads_to_assemble()  # generates an output

        # part2: create known output based on dummy data
        expected_output = 5

        # part3: testing the expected and test case
        self.assertEqual(test_my_function, expected_output)

    def test_2(self):
        """
        example of a negative test (exception, not equal)
        """
        # ---------------
        # exception way
        # ---------------

        # you are going to put arguments in here to cause the function to fail
        self.assertRaises(Exception, get_reads_to_assemble, arg1, arg2)

        # part1: your implementation output (test case)
        test_my_function = get_reads_to_assemble()  # generates an output

        # part3: testing the expected and test case
        self.assertEqual(test_my_function, expected_output)


"""
How to execute unit testing
- assumes that you have a test method


python -m unittest test_core.py
"""
