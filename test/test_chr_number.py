import sys

sys.path.append("../drive")
import utility_scripts

def test_chromosome_number():
    """Function to make sure that the proper chromosome number 
    is being returned"""
    
    test_str: str = "test_file_chr09.txt"

    chr_num: str = utility_scripts.get_chr_num(test_str, r'chr\d\d\.')

    assert chr_num == "chr09"

