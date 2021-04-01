import sys
sys.path.append("../drive")

from pre_shared_segments_analysis_scripts.shared_segment_detection.gathering_pairs.collect_shared_segments import generate_parameters, build_unique_id_dict, create_ibd_arrays

def test_generate_parameters():
    """unit test to test if the parameters are being properly generated"""

    # creating a list to keep track of errors
    errors: list = []

    param_dict: dict = generate_parameters("hapibd")

    # checking to make sure it returns the right type
    if type(param_dict) != dict:
        errors.append(f"Expected output to be of type dict, instead it was of type {type(param_dict)}")
    
    # making sure it is returning the right length. For 
    # hapibd it should be 6
    if len(param_dict) != 6:
        errors.append(f"Expected the resulting dictionary to have 6 key value pairs instead found {len(param_dict)}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_build_unique_id_dict():
    """unit test for the build_unique_id_dict"""

    # making an error log
    errors: list = []

    # making a list of iids
    iid_list: list = ["R123412", "R89767545", "R08912347"]

    # run the function to test
    uniq_dict: dict = build_unique_id_dict(iid_list)

    # checking how many keys the resulting dictionary has
    if len(uniq_dict) != 3:
        errors.append(f"Expected the resulting dictionary to have three keys, instead found {len(uniq_dict)}")
    # checking what the values of the dictionary is
    if [0,1,2] != list(uniq_dict.values()):
        errors.append(f"Expected the resulting dictionary to have the values 0,1, and 2 but instead found {list(uniq_dict.values())}")
    # checking the type of the dictionary
    if type(uniq_dict) != dict:
        errors.append(f"Expected the resulting dictionary to be of type dict, instead it was of type {type(uniq_dict)}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_create_ibd_arrays():
    """unit test to test the create_ibd_arrays"""

    # creating a list ot keep track of errors
    errors = []

    # running the create_ibd_arrays function
    IBD_data, IBD_index = create_ibd_arrays()

    # checking to make sure all the keys are right in the IBD_data dict
    if [str(i) for i in range(1,23)] != list(IBD_data.keys()):
        errors.append(f"Expected the dictionary to have keys from 1 to 23 but instead found {list(IBD_data.keys())}")
    # checking to make sure that a values in the dictionary are empty
    for i in IBD_data.values():
        if len(i) != 0:
            errors.append(f"Expected the values in the empty dictionary to be zero instead found a value with a length of {len(i)}")
    # checking to make sure all the keys are right in the IBD_index dict
    if [str(i) for i in range(1,23)] != list(IBD_index.keys()):
        errors.append(f"Expected the dictionary to have keys from 1 to 23 but instead found {list(IBD_index.keys())}")
    # checking to make sure the values in the IBD_index dictionary are right
    for i in IBD_index.keys():
        inner_dict: dict = IBD_index[i]
        if inner_dict["start"] != 999999999 and inner_dict["end"] != 0 and inner_dict["allpos"] != []:
            errors.append("The contents of the inner dictionary in the IBD_index dictionary are not as expected")   
      
    assert not errors, "errors occured: \n{}".format('\n'.join(errors))