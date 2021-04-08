import sys
import pandas as pd
import pytest
sys.path.append("../drive")

from pre_shared_segments_analysis_scripts.shared_segment_detection.gathering_pairs.find_IBD_phenotype import get_carriers, gather_gene_info, get_ibd_file

def test_get_carriers():
    """unit test for the get_carriers function"""
    # making a list to keep track of errors
    errors: list = []

    pheno_dict: dict = {
        "IID": ["R123541234", "R3456345", "R12039487", "R12341", "R12341234"]
        }

    pheno_df: pd.DataFrame = pd.DataFrame.from_dict(pheno_dict)

    iid_list: list = get_carriers(pheno_df)

    if len(iid_list) != 5:
        errors.append(f"Expected the output of the resulting iid_list to be 5 instead it was {len(iid_list)}")
    if type(iid_list) != list:
        errors.append(f"Expected the type of the iid_list to be a list instead it was {type(iid_list)}")

    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_get_carriers_len_exception():
    """unit test for checking the exception if the iid_list is empty in the get_carriers function"""
    
    pheno_dict: dict = {
        "IID": []
        }

    pheno_df: pd.DataFrame = pd.DataFrame.from_dict(pheno_dict)

    with pytest.raises(Exception) as e:
        assert get_carriers(pheno_df)
    assert str(e.value) == "Expected at least grid to be in the provided phenotype carriers file"

def test_get_carriers_keyerror():
    """unit test to check the keyerror exception if the dataframe does not have the IID column"""

    pheno_dict: dict = {
        "FID": ["R123541234", "R3456345", "R12039487", "R12341", "R12341234"]
        }

    pheno_df: pd.DataFrame = pd.DataFrame.from_dict(pheno_dict)    

    with pytest.raises(Exception):
        assert get_carriers(pheno_df)
    
def test_gather_gene_info():
    """unit test to test the gather_gene_info"""

    errors: list = []

    # running the function on the test data
    gene_dict: dict = gather_gene_info("./test_data/test_gmap_file.txt")

    # checking to see if the dictionary has the right keys
    if list(gene_dict.keys()) != ["MPO", "EXT1"]:
        errors.append(f"Expected the dictionary to have the keys 'MPO' and 'EXT1', instead the columns, {gene_dict.keys()} were found")
    # Next three lines check to make sure that the values for the inner dictionary from the M<P key are correct
    if int(gene_dict["MPO"]["chr"]) != 1:
        errors.append(f"Expected the value of the dictionary at ['MPO']['chr'] to be 1 instead it was {int(gene_dict['MPO']['chr'])}")
    
    if gene_dict["MPO"]["start"] != 1234:
        errors.append(f"Expected the value of the dictionary at ['MPO']['start'] to be 1234 instead it was {int(gene_dict['MPO']['start'])}")

    if gene_dict["MPO"]["end"] != 4567:
        errors.append(f"Expected the value of the dictionary at ['MPO']['end'] to be 4567 instead it was {int(gene_dict['MPO']['end'])}")
    
    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_get_ibd_file():
    """unit test for testing the get_ibd_file"""

    # creating a list to mimic a file list
    ibd_file_list: list = ["test_file_chr1.txt", "test_file_chr2.txt"]

    ibd_file: str = get_ibd_file(ibd_file_list, "1")
    
    assert ibd_file == "test_file_chr1.txt", f"Expected the file, test_file_chr1.txt, to be returned, instead the file was {ibd_file}"

def test_get_ibd_file_num_input():
    """unit test for testing the get_ibd_file when a integer is passed as chr_num"""

    # creating a list to mimic a file list
    ibd_file_list: list = ["test_file_chr1.txt", "test_file_chr2.txt"]

    ibd_file: str = get_ibd_file(ibd_file_list, 1)
    
    assert ibd_file == "test_file_chr1.txt", f"Expected the file, test_file_chr1.txt, to be returned, instead the file was {ibd_file}"