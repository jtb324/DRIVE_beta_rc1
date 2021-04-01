import sys
sys.path.append("../drive/")

import pre_shared_segments_analysis_scripts.generate_indx_dict as generate_indx_dict

def test_parameter_dict():
    """unit test for testing the Parameter_Dict class"""

    # creating a list to keep track of errors
    errors: list = []
    # generating the class object
    indx_dict_former = generate_indx_dict.Parameter_Dict("hapibd")

    indx_dict: dict = indx_dict_former.return_param_dict()

    # checking to make sure the ibd_program is correct
    if indx_dict_former.ibd_program != "hapibd":
        errors.append(f"Expected the ibd_program attribute to be hapibd. Instead it was {indx_dict_former.ibd_program}")
    # checking to make sure the parameters dictionary has the right number of keys
    if len(indx_dict_former.param_dict) != 5:
        errors.append(f"Expected the dictionary to have 5 keys in the base class. Instead found {len(indx_dict_former.param_dict)}")
    # checking to make sure the function returns a dictionary as 
    # expected
    if type(indx_dict) != dict:
        errors.append(f"Expected the return_param_dict method to return a dictionary, instead returned {type(indx_dict)}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_germline_indices():
    """unit test for the Germline_Indices class"""
    # creating a list to keep track of errors
    errors: list = []
    # generating the class object
    indx_dict_former = generate_indx_dict.Germline_Indices("germline")

    indx_dict_former.update_indices()

    indx_dict: dict = indx_dict_former.return_param_dict()

    # checking to make sure the ibd_program is correct
    if indx_dict_former.ibd_program != "germline":
        errors.append(f"Expected the ibd_program attribute to be germline. Instead it was {indx_dict_former.ibd_program}")
    # checking to make sure the parameters dictionary has the right number of keys
    if len(indx_dict_former.param_dict) != 7:
        errors.append(f"Expected the dictionary to have 7 keys in the base class. Instead found {len(indx_dict_former.param_dict)}")
    
    # checking to make sure the function returns a dictionary as 
    # expected
    if type(indx_dict) != dict:
        errors.append(f"Expected the return_param_dict method to return a dictionary, instead returned {type(indx_dict)}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_hapibd_indices():
    """unit test for the Hapibd_Indices class"""
    # creating a list to keep track of errors
    errors: list = []
    # generating the class object
    indx_dict_former = generate_indx_dict.Hapibd_Indices("hapibd")

    indx_dict_former.update_indices()

    indx_dict: dict = indx_dict_former.return_param_dict()

    # checking to make sure the ibd_program is correct
    if indx_dict_former.ibd_program != "hapibd":
        errors.append(f"Expected the ibd_program attribute to be hapibd. Instead it was {indx_dict_former.ibd_program}")
    # checking to make sure the parameters dictionary has the right number of keys
    if len(indx_dict_former.param_dict) != 6:
        errors.append(f"Expected the dictionary to have 6 keys in the base class. Instead found {len(indx_dict_former.param_dict)}")
    
    # checking to make sure the function returns a dictionary as 
    # expected
    if type(indx_dict) != dict:
        errors.append(f"Expected the return_param_dict method to return a dictionary, instead returned {type(indx_dict)}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_ilash_indices():
    """unit test for the Ilash_Indices class"""
    # creating a list to keep track of errors
    errors: list = []
    # generating the class object
    indx_dict_former = generate_indx_dict.Ilash_Indices("ilash")

    indx_dict_former.update_indices()

    indx_dict: dict = indx_dict_former.return_param_dict()

    # checking to make sure the ibd_program is correct
    if indx_dict_former.ibd_program != "ilash":
        errors.append(f"Expected the ibd_program attribute to be ilash. Instead it was {indx_dict_former.ibd_program}")
    # checking to make sure the parameters dictionary has the right number of keys
    if len(indx_dict_former.param_dict) != 6:
        errors.append(f"Expected the dictionary to have 6 keys in the base class. Instead found {len(indx_dict_former.param_dict)}")
    
    # checking to make sure the function returns a dictionary as 
    # expected
    if type(indx_dict) != dict:
        errors.append(f"Expected the return_param_dict method to return a dictionary, instead returned {type(indx_dict)}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))
