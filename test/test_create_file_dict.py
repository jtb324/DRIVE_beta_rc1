import sys
import pandas as pd
import collections
sys.path.append("../drive")

from pre_shared_segments_analysis_scripts.file_dict_creator.create_file_dict import form_chr_num, find_match, get_file_dict

def test_form_chr_num():
    """unit test for the form_chr_num function which returns the proper 
    chromosome number to be used as a key in a dictionary"""
    # creating a list to keep all the errors
    errors: list = []

    chr_num_1: str = form_chr_num(9)
    chr_num_2: str = form_chr_num(11)

    if chr_num_1 != "chr09":
        errors.append(f"The expected output of the form_chr_num function was chr09, instead the function returned {chr_num_1}")

    if chr_num_2 != "chr11":
        errors.append(f"The expected output of the form_chr_num function was chr10, instead the function returned {chr_num_2}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_find_match_carrier():
    """unit test to test the find_match function"""

    file_list: list = [
        "test_file_1_chr01.single_var.csv",
        "test_file_2_chr10.single_var.csv",
        "test_file_3_chr19.single_var.csv"
    ]

    assert "chr19" in find_match(file_list, "chr19", "carrier")

def test_find_match_ibd():
    """unit test to test the find_match function"""
    errors: list = []

    file_list: list = [
        "test_file_1_chr1.match.gz",
        "test_file_2_chr10.match.gz",
        "test_file_3_chr19.match.gz"
    ]

    if not "chr1." in find_match(file_list, "chr01", "ibd"):
        errors.append("chr1. was not found within the output file")

    if not "chr10" in find_match(file_list, "chr10", "ibd"):
        errors.append("chr10 was not found with the output file string")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))


def test_find_match_map():
    """unit test to test the find_match function"""

    file_list: list = [
        "test_file_1_chr01_variants.map",
        "test_file_2_chr10_variants.map",
        "test_file_3_chr19_variants.map"
    ]

    assert "chr10" in find_match(file_list, "chr10", "map")

def test_get_file_dict():
    """unit test to test the get_file_dict function"""
    # making a log to keep the errors
    errors: list = []

    # making three sets of list for the test
    map_file_list: list = ["test_file_1_chr01_variants.map"]
    ibd_file_list: list = ["test_file_1_chr1.match.gz"]
    carrier_file_list: list = ["test_file_1_chr01.single_var.csv"]

    file_dict: dict = get_file_dict(
        map_file_list, carrier_file_list, ibd_file_list)
    # if statement to check if the inner dictionary is the correct size
    if len(file_dict["chr01"]) != 3:
        errors.append(
            "The dictionary formed has the expected n7umber of files for chr01")
    # checking to make sure that the string for hte value
    #  connected to the carrier key is the right string
    if file_dict["chr01"]["carrier"] != carrier_file_list[0]:
        errors.append("The expected file {} was found within the corresponding key in the output dictionary".format(
            carrier_file_list[0]))
    # Making sure that all the keys are present in the output 
    # dictionary
    if collections.Counter(list(file_dict["chr01"].keys())) != collections.Counter(["carrier", "map", "ibd"]):
        errors.append("expected the output dictionary to have three keys carrier, ibd, and map. Found keys: {} instead.".format(list(file_dict["chr01"].keys())))

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))