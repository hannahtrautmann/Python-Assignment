from python_assignment import *
import pytest

def test_observed_kmer_dict_single_value():
    """
    Summary line: tests observed_kmer_dict_single_value function

    Extended description: uses an example sequence and an example k to determine if appropriate kmer dictionary is made
    """
    example_sequence = 'ATTTGGATT'
    example_k = 3
    expected_dict = {'ATT': 2, 'TTT': 1, 'TTG': 1, 'TGG': 1, 'GGA': 1, 'GAT': 1}
    assert expected_dict == observed_kmer_dict_single_value(example_sequence, example_k)

def test_observed_kmers_list():
    """
    Summary line: tests observed_kmer_list function

    Extended description: uses an example sequence to determine if appropriate kmer dictioary is made for each possible k,
    with correct number of values for observed kmers.
    """
    example_sequence = 'ATTTGGATT'
    expected_list = [3,5,6,6,5,4,3,2,1]
    assert expected_list == observed_kmers_list(example_sequence)

def test_observed_kmers_list_empty():
    """
    Summary line: tests observed_kmer_list function with empty sequence

    Extended description: uses an empty sequence to determine if function will return an empty list.
    """
    example_sequence = ''
    expected_list = []
    assert expected_list == observed_kmers_list(example_sequence)

def test_possible_kmers_list():
    """
    Summary line: tests possible_kmer_list function

    Extended description: uses an example sequence to determine if appropriate kmer dictioary is made for each possible k,
    with correct number of values for possible kmers.
    """
    example_sequence = 'ATTTGGATT'
    expected_list = [4,8,7,6,5,4,3,2,1]
    assert expected_list == possible_kmers_list(example_sequence)

def test_possible_kmers_list_empty():
    """
    Summary line: tests possible_kmer_list function with empty sequence

    Extended description: uses an empty sequence to determine if function will return an empty list.
    """
    example_sequence = ''
    expected_list = []
    assert expected_list == possible_kmers_list(example_sequence)

def test_kmers_df():
    """
    Summary line: tests kmers_df function

    Extended description: uses an example sequence to determine if correct dataframe is made for k values, observed kmers, and possible kmers.
    """
    example_sequence = 'ATTTGGATT'
    expected_df = pd.DataFrame(
        {'k':[1,2,3,4,5,6,7,8,9], 'observed_kmers':[3,5,6,6,5,4,3,2,1], 'possible_kmers':[4,8,7,6,5,4,3,2,1]},\
        index = [0,1,2,3,4,5,6,7,8]
        )
    assert expected_df.equals(kmers_df(example_sequence))

def test_empty_kmers_df():
    """
    Summary line: tests kmers_df function with empty sequence

    Extended description: uses an empty sequence to determine if function will return an empty database.
     """
    example_sequence = ''
    expected_df = pd.DataFrame(
        {'k':[], 'observed_kmers':[], 'possible_kmers':[]},\
        index = []
        )
    assert expected_df.equals(kmers_df(example_sequence))

def test_linguistic_complexity():
    """
    Summary line: tests linguistic_complexity function with example sequence

    Extended description: uses an example sequence to determine if function will accurately calculate linguistic complexity.
     """
    example_sequence = 'ATTTGGATT'
    expected_output = 0.875
    assert expected_output == linguistic_complexity(example_sequence)
