import sys
sys.path.append("../drive")

from pre_shared_segments_analysis_scripts.Shared_segment_Formatter_CLI import find_match

def test_find_match_carrier():
    """unit test to test the find_match function"""

    file_list: list = [
        "test_file_1_chr01.single_var.csv", 
        "test_file_2_chr10.single_var.csv",
        "test_file_3_chr19.single_var.csv"
        ]
    
    assert "chr19" in find_match(file_list, "chr19", "carrier")

def test_find_match_ibd():
    """unit test to test the find_match function"""
    errors: list = []

    file_list: list = [
        "test_file_1_chr1.match.gz", 
        "test_file_2_chr10.match.gz",
        "test_file_3_chr19.match.gz"
        ]

    if not "chr1." in find_match(file_list, "chr01", "ibd"):
        errors.append("chr1. was not found within the output file")
    
    if not "chr10" in find_match(file_list, "chr10", "ibd"):
        errors.append("chr10 was not found with the output file string")
    
    assert not errors, "errors occured: \n{}".format('\n'.join(errors))

def test_find_match_map():
    """unit test to test the find_match function"""

    file_list: list = [
        "test_file_1_chr01_variants.map", 
        "test_file_2_chr10_variants.map",
        "test_file_3_chr19_variants.map"
        ]
    
    assert "chr10" in find_match(file_list, "chr10", "map")
