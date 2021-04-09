import sys
import os
import pandas as pd
import collections
sys.path.append("../drive")

from pre_shared_segments_analysis_scripts.shared_segment_detection.gathering_pairs.gather_ibd_info import create_iid_dict, create_dict_with_var_pos, create_no_carriers_file, create_var_info_dict, filter_no_carriers


def test_create_iid_dict():
    """unit test for the create_iid_dict to make sure it is adding the value to the dictionary"""

    # making a log to keep the errors
    errors: list = []

    variant: str = "exm65234"

    iid_dict: dict = {}

    carrier_dict: dict = {
        "IID": ["R123541334", "R098764567", "r1234568"],
        "Variant ID":["exm65234", "exm65234", "exm65234"]
    }

    carrier_df: pd.DataFrame = pd.DataFrame.from_dict(carrier_dict)

    iid_dict = create_iid_dict([variant], iid_dict, carrier_df)
    print(iid_dict)

    if len(iid_dict[variant]) != 3:
        errors.append(f"expected the length of the value for the key {variant} to be 3, instead the length was {len(iid_dict[variant])}")
    
    if variant not in iid_dict.keys():
        errors.append(f"Variant, {variant}, was expected to be in the key list of the iid_dict. The value was not found")
    
    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_create_dict_with_var_pos():
    """unit test for the create_dict_with_var_pos function to make sure it is also adding the value correctly"""
    # making a log to keep the errors
    errors: list = []

    # Getting the necessary parameters for the function
    variant: str = "exm65234-G"

    var_pos_dict: dict = {}

    map_dict: dict = {
        "chr":[10, 10],
        "variant id":["exm65234", "var123453124"],
        "cM": [0,0],
        "site": [12341234, 34563456]
    }

    map_df: pd.DataFrame = pd.DataFrame.from_dict(map_dict)

    # running the create_dict_with_var_pos functions
    var_pos_dict = create_dict_with_var_pos([variant], var_pos_dict, map_df)

    if type(var_pos_dict[variant]) == str:
        errors.append(f"expected the length of the value for the key {variant} to be 3, instead the length was {len(var_pos_dict[variant])}")
    
    if variant not in var_pos_dict.keys():
        errors.append(f"Variant, {variant}, was expected to be in the key list of the iid_dict. The value was not found")
    
    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_create_no_carriers_file():
    """unit test to check if the file no_carriers_in_file.txt is being formed"""

    create_no_carriers_file("exm3425234", "chr12", "./")

    assert os.path.exists("./no_carriers_in_file.txt"), "file was not formed properly"

    if os.path.exists("./no_carriers_in_file.txt"):
        os.remove("./no_carriers_in_file.txt")

def test_create_var_info_dict():
    """unit test for testing the variant info dict"""
    # making a log to keep the errors
    errors: list = []
    # creating the necessary parameters for the test
    var_info_dict: dict = {}

    var_iid_dict: dict = {
        "exm6523": ["R1234123", "R457645"]
    }
    var_dict: dict = create_var_info_dict(var_info_dict, var_iid_dict, "exm6523", "62345234")

    if "base_pos" not in var_dict["exm6523"].keys() or "iid_list"  not in var_dict["exm6523"].keys():
        errors.append("Expected to find two keys, base_pos and iid_list in the output keys. These keys were not found")
    if "exm6523" not in var_dict:
        errors.append("Expected to find the key exm6523 in the resulting dictionary but did not find this key")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_filter_no_carriers():
    """unit test for testing the filter_no_carriers function"""

    # making a log to keep the errors
    errors: list = []

    # creating a list of dictionary of variants and iids that carry the variants
    var_iid_dict: dict = {
        "exm6523": ["R1234123", "R457645"],
        "var09i21": ["R1234123", "R0897"],
        "10:1234124":["None"]
    }
    # running the function
    filter_dict: dict = filter_no_carriers(var_iid_dict, "./", "chr09")

    if len(filter_dict.keys()) != 2:
        errors.append(f"Expected the resulting dictionary to have two keys. {len(filter_dict.keys())} were found instead")
    
    if not os.path.exists("./no_carriers_in_file.txt"):
        errors.append("./no_carriers_in_file.txt, file was not formed properly")
    else:
        os.remove("./no_carriers_in_file.txt")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))


