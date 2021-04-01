import sys
sys.path.append("../drive")

from pre_shared_segments_analysis_scripts.IBDinput_byBP_classes import generate_parameters

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