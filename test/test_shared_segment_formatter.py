import sys
import os
import pandas as pd
import collections
sys.path.append("../drive")

from pre_shared_segments_analysis_scripts.Shared_segment_Formatter_CLI import create_iid_dict, create_dict_with_var_pos, create_no_carriers_file


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

    iid_dict = create_iid_dict(variant, iid_dict, carrier_df)

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
    variant: str = "exm65234"

    var_pos_dict: dict = {}

    map_dict: dict = {
        "chr":[10, 10],
        "variant id":["exm65234", "var123453124"],
        "cM": [0,0],
        "site": [12341234, 34563456]
    }

    map_df: pd.DataFrame = pd.DataFrame.from_dict(map_dict)

    # running the create_dict_with_var_pos functions
    var_pos_dict = create_dict_with_var_pos(variant, var_pos_dict, map_df)

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
