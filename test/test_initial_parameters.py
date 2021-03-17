# Setting the system path
import logging
import sys
sys.path.append('../drive')

from utility_scripts.user_input.initial_parameters import Input_Gather, get_dict_of_variables
# Importing modules

input_gather_object = Input_Gather()


def test_class_type():
    """Function to make suerthe class is of the expected type"""

    assert type(input_gather_object) == Input_Gather


def test_analysis_type():
    """Function to make sure that the analysis type will be caught if it is not a recognized argument"""

    assert input_gather_object.ANALYSIS_TYPE in ["gene", "maf", ""]


def test_min_cm():
    """This value test to make sure that the minimum centimorgan value is a integer"""

    assert type(input_gather_object.MIN_CM) == int


def test_threads():
    """Function checks to make sure that the number of threads passed to the program is an integer"""

    assert type(input_gather_object.THREADS) == int


def test_maf_filter():
    """unit test to make sure that the minor allele frequency threshold is a string"""
    assert type(input_gather_object.MAF_THRESHOLD) == str


def test_dictionary_keys():
    """unit test to make sure that all the expected dictionary keys are in the returned dictionary"""
    parameters_dict: dict = get_dict_of_variables(input_gather_object)

    para_keys: list = list(parameters_dict.keys())

    assert "ANALYSIS_TYPE" in para_keys and "MIN_CM" in para_keys and "THREADS" in para_keys and "MAF_THRESHOLD" in para_keys
