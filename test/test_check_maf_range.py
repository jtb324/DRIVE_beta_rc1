import sys
import pandas as pd
import pytest
sys.path.append("../drive")

from carrier_analysis_scripts.check_maf_range import load_frequency_file, filter_for_higher_maf, get_filtered_var_list, check_mafs

def test_load_frequency_file_success():
    """unit test for load_frequency_file"""
    freq_df: pd.DataFrame = load_frequency_file("./test_data/test_allele_freq.txt")

    assert type(freq_df) == pd.DataFrame, "The resulting output is not of the expected type pd.DataFrame instead it is: {1}".format(type(freq_df))

def test_load_frequency_file_exception():
    """unit test for the load_frequency_file function when an exception is raised"""
    with pytest.raises(Exception) as e:
        assert load_frequency_file("file.txt")
    assert str(e.value) == "file was not found"
    

    