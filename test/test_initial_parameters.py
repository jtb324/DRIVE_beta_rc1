# Setting the system path
import logging
import sys
from typing import List
sys.path.append('../drive')

from utility_scripts.user_input.initial_parameters import  check_provided_maf, check_provided_integer
# Importing modules

def test_check_provided_maf():
    """unit test to check if the returned maf is correct"""
    error_list: List[str] = []
    maf_threshold: float = check_provided_maf(0.10)

    if type(maf_threshold) != float:
        error_list.append(f"The resulting output is not of the expected type float instead it is: {(type(maf_threshold))}")
    if maf_threshold != 0.10:
        error_list.append(f"Expected the output to be 0.10 instead it was {maf_threshold}")

    assert not error_list, "errors occured: \n{}".format('\n'.join(error_list))

def test_check_provided_integer():
    """unit test to check if the returned maf is correct"""
    error_list: List[str] = []
    return_integer: int = check_provided_integer(3)

    if type(return_integer) != int:
        error_list.append(f"The resulting output is not of the expected type float instead it is: {(type(return_integer))}")

    if return_integer != 3:
        error_list.append(f"Expected the output to be 0.10 instead it was {return_integer}")

    assert not error_list, "errors occured: \n{}".format('\n'.join(error_list))