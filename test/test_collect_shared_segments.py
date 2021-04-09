import sys
import pandas as pd
import numpy as np
sys.path.append("../drive")

from pre_shared_segments_analysis_scripts.shared_segment_detection.gathering_pairs.collect_shared_segments import generate_parameters, build_unique_id_dict, create_ibd_arrays, get_pair_string, build_ibddata_and_ibddict, filter_for_gene_site

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

def test_get_pair_string():
    """unit test for testing the get_pair_string function"""

    # making a list to keep track of all the errors
    errors: list = []
    # testing the case where both ids are in the uniqID dict and they are in order
    row: pd.Series = pd.Series(np.array(["R12341234", "R2345112", "3.4"]))

    uniqID: dict = {
        "R12341234": 0,
        "R2345112": 1
    }

    both_pair_str: str = get_pair_string(row, 0,1,2, uniqID)

    if both_pair_str != '{0}:{1}-{2}'.format("3.4", "R12341234", "R2345112"):
        errors.append(f"Expected the string to be of the format: 3.4: R12341234-R2345112, instead it was of the format: {both_pair_str}")

    # testing the case where both ids are in the uniqID
    uniqID_2:dict = {
        "R12341234": 1,
        "R2345112": 0
    }

    out_of_order_str: str = get_pair_string(row, 0, 1, 2, uniqID_2)

    if out_of_order_str != '{0}:{1}-{2}'.format("3.4",  "R2345112", "R12341234"):
        errors.append(f"Expected the string to be of the format: 3.4: R2345112-R12341234, instead it was of the format: {out_of_order_str}")

    # testing the case where 1st id is in uniqID. Use uniqID 
    row1: pd.Series = pd.Series(np.array(["R12341234", "R234545678", "3.4"]))

    first_in_pair: str = get_pair_string(row1, 0,1,2, uniqID)

    if first_in_pair != '{0}:{1}-{2}'.format("3.4",  "R12341234", "R234545678"):
        errors.append(f"Expected the string to be of the format: 3.4: R12341234-R234545678, instead it was of the format: {first_in_pair}")

    # testing the case where 2 id is in the uniqID. Use uniqID
    row2: pd.Series = pd.Series(np.array(["R124765234", "R2345112", "3.4"]))
    
    second_in_pair: str = get_pair_string(row2, 0, 1, 2, uniqID)

    if second_in_pair != '{0}:{1}-{2}'.format("3.4", "R2345112", "R124765234"):
        errors.append(f"Expected the string to be of the format: 3.4: R2345112-R124765234, instead it was of the format: {second_in_pair}")
    
    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_apply_on_get_pair_string():
    """unit test to make sure that the apply function will work on the get_pair_str """
    # making a list to keep track of errors
    errors: list = []

    uniqID: dict = {
        "R12341234": 0,
        "R2345112": 2,
        "R1234512341":1,
        "R1290384": 3
    }

    chunk_dict: dict = {
        0: ["R12341234", "R2345112", "R1234512341", "4356789"],
        1: ["R2345112", "R1234512341", "R097897", "R1290384"], 
        2: ["4.5", "5.5", "6.6", "7.7"]
    }
    chunk: pd.DataFrame = pd.DataFrame.from_dict(chunk_dict)

    chunk["pair_string"] = chunk.apply(lambda row: get_pair_string(row, 0, 1, 2, uniqID), axis = 1)

    # checking if the pair_string column is formed
    if "pair_string" not in chunk.columns:
        errors.append(f"expected the dataframe to have a column titled 'pair_string', but this was not found. Instead the columns are {chunk.columns}")

    # checking if a value in the column is as expected
    if chunk[chunk[0] == "R12341234"]["pair_string"].values.tolist()[0] != '4.5:R12341234-R2345112':
        errors.append(f"Expected the pair string to be: '4.5:R12341234-R2345112', instead it was {chunk[chunk[0] == 'R12341234']['pair_string'].values.tolist()[0]}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_build_ibddata_and_ibdindex():
    """unit test to test the build_ibddata_and_ibdindex function"""

    # making a list to keep track of all the errors
    errors: list = []
    # testing the case where both ids are in the uniqID dict and they are in order
    row: pd.Series = pd.Series(np.array(["R12341234", "R2345112","4", "3.4","12341234", "23452345",]))
    
    row = row.append(pd.Series(["3.4:R12341234-R2345112"], index=["pair_string"]))

    IBDdata: dict = {str(i): {} for i in range(1, 23)}
    IBDindex: dict = {
            str(i): {
                'start': 999999999,
                'end': 0,
                'allpos': []
            }
            for i in range(1, 23)
        }
    
    chr_num: pd.Series = build_ibddata_and_ibddict(row, 4,5,2, IBDdata, IBDindex)

    # checking if a value got put into the IBDdata
    if len(IBDdata["4"]) == 0:
        errors.append(f"Expected the IBDdata key 4 to not be empty.")
    if len(IBDindex["4"]["allpos"]) == 0:
        errors.append(f"Expected the allpos position for chromosome 4 to not be empty")
    if chr_num != "4":
        errors.append(f"Expected chromosome number returned was 4, instead returned {chr_num}")
    
    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_filter_for_gene_site():
    """Unit test for the filter_for_gene_site function"""

    # creating a list to keep track of errors
    errors: list = []

    gene_start: int = 20
    gene_end : int = 50

    chunk_dict: dict = {
        0:["R234523", "R45674567", "R234523444", "R23452345", "R98766"],
        1:["R45674567", "R234523444", "R23452345", "R98766", "R89"],
        2: ["12", "23", "10", "55", "10"],
        3:["40", "55", "56", "70", "15"]
    }    

    chunk_df: pd.DataFrame = pd.DataFrame.from_dict(chunk_dict)

    filtered_df: pd.DataFrame = filter_for_gene_site(chunk_df, 2,3, gene_start, gene_end)
    
    if len(filtered_df) != 3:
        errors.append(f"Expected the length of the returned dataframe to be 3 instead it was {len(filtered_df)}")
    if (filtered_df[0].values != ["R234523", "R45674567", "R2345234"]).all():
        errors.append(f"Expected the returned dataframe to contain the pair 1 values of ['R234523', 'R45674567', 'R234523'], instead the values {filtered_df[0].values}")
    if (filtered_df[2].values != ["12", "23", "10"]).all():
        errors.append(f"Expected the returned dataframe to have the values ['12', '23', '10'] for the start_indx column instead found the values {filtered_df[2].values}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))