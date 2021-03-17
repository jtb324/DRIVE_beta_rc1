
import sys
sys.path.append("../drive")

from carrier_analysis_scripts.identify_single_var_carrier import find_all_files
import file_creator_scripts
import population_filter_scripts
import utility_script

# secondary imports


def test_all_files_same_extension():
    """This function will make sure that the output of 
    the find_all_files function has the same extensions"""

    output_file_list: list = find_all_files("./test_data")

    checked_list: list = [
        file for file in output_file_list if file[-4:] == ".raw"]

    assert len(checked_list) == len(output_file_list)
