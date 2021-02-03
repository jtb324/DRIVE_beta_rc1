import analysis_haplotypes
import os
import sys
from Levenshtein import distance as levenshtein_distance


def check_if_equal(pair_1_str: str, pair_2_str: str) -> int:
    '''This function will check to see if the strings are equal'''

    if pair_1_str == pair_2_str:
        return 1
    else:
        return 0


def compare_haplotype_str(file_tuple: tuple) -> int:
    '''This function will compare the haplotype strings and then write that to file'''

    pair_1_haplotype_str: str = file_tuple[0].split(" ", 6)[-1]

    pair_2_haplotype_str: str = file_tuple[1].split(" ", 6)[-1]

    # using a unit test to make sure the strings are the same length
    analysis_haplotypes.check_strs_len(pair_1_haplotype_str,
                                       pair_2_haplotype_str)
    print(len(pair_1_haplotype_str))

    # check to see if the strings are the same. If they are the value will be 1
    are_equal_int: int = check_if_equal(pair_1_haplotype_str,
                                        pair_2_haplotype_str)

    if are_equal_int:

        haplotype_diff: int = 0

    else:
        # These next few lines initialize the levenshtein distance array

        haplotype_diff = levenshtein_distance(
            pair_1_haplotype_str, pair_2_haplotype_str)

    return haplotype_diff
