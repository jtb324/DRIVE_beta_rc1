import sys
import pandas as pd
sys.path.append("../drive")

from pre_shared_segments_analysis_scripts.shared_segment_detection.gathering_pairs.filtering_functions import filter_to_individual_in_uniqID, filter_to_greater_than_3_cm, filter_for_correct_base_pair

def test_filter_to_individual_in_uniqID():
    """unit test to test the filter_to_individual_in_uniqID function"""

    chunk_info: dict = {
        0: ["ind1", "ind2"],
        1: ["pair1", "pair2"],
    } 
    chunk: pd.DataFrame = pd.DataFrame.from_dict(chunk_info)

    uniqID: dict = {
       "ind1":0,
       "ind3":1     
    }

    chunk_subset: pd.DataFrame = filter_to_individual_in_uniqID(chunk, uniqID,0, 1)

    assert len(chunk_subset) <= len(chunk), "Expected the length of the resulting dataframe to be less than or equal to the length of the original."

def test_filter_to_greater_than_3_cm():
    """unit test to test the filter_to_greater_than_3_cm"""
    # creating a list to keep track of errors
    errors = []

    # creating the dataframe chunk that is needed
    chunk_info: dict = {
        0: ["ind1", "ind2"],
        1: ["pair1", "pair2"],
        2: [2.91, 3.9]
    } 

    chunk: pd.DataFrame = pd.DataFrame.from_dict(chunk_info)

    chunk_subset: pd.DataFrame = filter_to_greater_than_3_cm(chunk, 2, 3)

    # checking to make sure that the resulting dataframe is not larger than the initial dataframe
    if len(chunk_subset) > len(chunk):
        errors.append(f"Expected maximum length of the resulting dataframe: {len(chunk)} \n Length of resulting dataframe was instead {len(chunk_subset)}")
    # checking to make sure that it found the right row
    if chunk_subset[0].values.tolist()[0].strip() != "ind2":
        errors.append(f"Expected the function to return the row for ind2, instead returned the row for {chunk_subset[0].values.tolist()[0]}")
    # checking to see if the value in column 2 of the chunk subset is greater than the threshold of 3 in this case
    if chunk_subset[2].values.tolist()[0] < 3:
        errors.append(f"Expected the min_CM threshold to be greater than 3 instead it was {chunk_subset[2]}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_filter_for_correct_base_pair():
    """unit test for testing the filter_for_correct_base_pair function"""
    # creating a list to keep track of the errors
    errors = []

    # creating the dataframe chunk that is needed
    chunk_info: dict = {
        0: ["ind1", "ind2", "ind3"],
        1: ["pair1", "pair2", "pair3"],
        2: [2.91, 3.9, 4],
        3: [5,8, 5],
        4: [10,11, 6]
    } 

    chunk: pd.DataFrame = pd.DataFrame.from_dict(chunk_info)

    chunk_subset: pd.DataFrame = filter_for_correct_base_pair(chunk, 3, 4, 7)

    # checking to make sure that the resulting dataframe is not larger than the initial dataframe
    if len(chunk_subset) > len(chunk):
        errors.append(f"Expected maximum length of the resulting dataframe: {len(chunk)} \n Length of resulting dataframe was instead {len(chunk_subset)}")
    # checking to make sure that it found the right row
    if chunk_subset[0].values.tolist()[0].strip() != "ind1":
        errors.append(f"Expected the function to return the row for ind2, instead returned the row for {chunk_subset[0].values.tolist()[0]}")
    # checking to see if the value in column 3 is less than 7
    if not chunk_subset[3].values.tolist()[0] < 7:
        errors.append(f"Expected the start position to be less than seven instead it was {chunk_subset[3].values.tolist()[0]}")
    # checking to see if the value in column 2 of the chunk subset is greater than the threshold of 3 in this case
    if not chunk_subset[4].values.tolist()[0] > 7:
        errors.append(f"Expected the end position to be greater than seven instead it was {chunk_subset[4].values.tolist()[0]}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))