# This file wil contain unit test that test certain aspects of the files in this directory


def check_tuple_size(file_tuple: tuple):
    '''This unit test will check the size of the tuples'''

    assert (len(file_tuple)
            ) == 2, "The tuple is unexpectedly longer than two elements"


def check_strs_len(str_1: str, str_2: str):
    '''This unit test will check to make sure the two strings for each pair are the same len'''
    assert (
        len(str_1) == len(str_2)
    ), "The strings for both individuals in the pair are not the same length"
