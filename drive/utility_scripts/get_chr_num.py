import re


def get_chr_num(file: str, pattern: str) -> str:
    """This function will get the chr_num from the file name
    Parameters
    __________
    file : str
        string containing the name of the file that the chromosome 
        number will be extracted from
        
    pattern : str
        this is the regular expression that will be matched
        
    Returns
    _______
    str
        returns a string of the chromosome number. Should be of 
        the format chrXX where X is a number 
    """

    match = re.search(pattern, file)

    chr_num = match.group(0)

    chr_num = chr_num.strip(".").strip("_")

    return chr_num


def get_alt_chr_num(chr_num: str) -> str:
    '''This function is used to convert the chromosome numbers 1-9 from
    a format of .chr01_ to .chr1_ incase the file needs that. Some older
    files may use this format and this will allow the program to keep
    running'''

    return "".join(chr_num.split("0"))


def add_zero_to_chr_num(chr_num: str) -> str:
    '''This function will take the chromsome number format of chr1 and convert to 
    chr01 for all single digit numbers'''

    chr_str: str = chr_num[0:3]

    chr_number: str = chr_num[-1:]

    return "".join([chr_str, "0", chr_number])

def match_chr(search_list: list, string: str) -> str:
    """Function that will find the chromosome number of the string
    Parameters
    __________
    search_list : list
        list of the different regular expressions to search for
    
    string : str
        string to match a pattern against
    Returns
    ________
    str
        string containing the chromosome number that was matched"""
    for search_term in search_list:

        match = re.search(search_term, string)

        if match:
            
            return match.group(0)
