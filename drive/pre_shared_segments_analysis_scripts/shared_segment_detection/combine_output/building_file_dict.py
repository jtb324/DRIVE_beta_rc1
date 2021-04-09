# this script helps build the dictionary that contains keys of unique information like the chromosome number and the variant id/gene name (depends on the analysis type) and then the values are the .small.txt.gz files that are lined up with the correct ibd program

import re

def get_variant_id(ibd_file: str) -> str:
    """This function will get the proper variant_id
    Parameters
    __________
    ibd_file : str
        This string list the ibd_file of interest
    
    Returns 
    _______
    str
        returns the variant_id as a string"""

    dot_indx_list: list = [
        indx for indx, letter in enumerate(ibd_file) if letter == "."
    ]

    # this gets the string up until the dot_indx_list value by 
    variant_id_handler = {
        4: ibd_file[:dot_indx_list[0]],
        5: ibd_file[:dot_indx_list[1]]
    }

    variant_id: str = variant_id_handler[len(dot_indx_list)]

    return variant_id

def get_gene_name(file_name: str) -> str:
    """Function to get the gene name from the file name
    Parameters
    __________
    file_name : str
        This is the file name of the .small.txt.gz file
    
    Returns
    _______
    str
        returns a string of just the gene name
    """
    # getting a list of all dot indices
    dot_indx_list: list = [
        indx for indx, letter in enumerate(file_name) if letter == "."
    ]

    # getting the gene name by getting the substring that is up until the first .
    return file_name[:dot_indx_list[0]]

def match_chr(search_list: list, file: str) -> str:
    """Function that will find the chromosome number of the string
    Parameters
    __________
    search_list : list
        list of the different regular expressions to search for
    
    file : str
        string to match a pattern against
    Returns
    ________
    str
        string containing the chromosome number that was matched"""
    for search_term in search_list:

        match = re.search(search_term, file)

        if match:
            
            return match.group(0)

def get_ibd_file_substring(program_list: list, ibd_file: str) -> str:
    """Function to find the substring of the ibd_program from the ibd_file
    Parameters
    __________
    program_list : list
        the is a list of the ibd programs and will be either ilash 
        or hapibd
    
    ibd_file : str
        this is the string of the .small.txt.gz files
    
    Returns
    _______
    str
        returns a string that list the shorted ibd file from after 
        the words of either hapibd or ilash
    """
    # iterating through the program list which is either ilash or hapibd
    for ibd_program in program_list:

        if ibd_program in ibd_file:
            underscore_indx = ibd_file.find("".join([ibd_program, "_"]))

            shorten_ibd_file_string: str = ibd_file[underscore_indx + len(ibd_program) + 1:]
    
    return shorten_ibd_file_string

def determine_identifier(short_ibd_string: str, analysis_type) -> str:
    """Function to return either the variant id or the gene 
    depending on what analysis mode is provided
    Parameters
    __________
    short_ibd_string : str
        string that has the .small.txt.gz file without the ibd_progm 
        name
    analysis_type : str
        this is the provide analysis type string
    
    Returns
    _______
    str
        returns a string of either the variant it or the gene name
    """

    identifier_handler: dict = {
        True: get_gene_name(short_ibd_string),
        False: get_variant_id(short_ibd_string)
    }

    return identifier_handler[analysis_type=="phenotype"]

def build_file_dict(ibd_file_list: list, program_list: list, analysis_type: str) -> dict:
    """This function gathers the necessary files per chromosome and variant/gene combo
    Parameters
    __________
    ibd_file_list : list
        This is a list of all the .small.txt.gz files
        
    program_list : list
        This is a list of ibd detection programs used. Right now 
        it will be either hapibd or ilash
        
    analysis_type : str
        string that tells what type of analysis is being used. If 
        the phenotype analysis is used the code will branch one way 
        else it will go a different way
    
    Returns
    _______
    dict
        returns a dictionary where the keys are either a tuple of chr_num and variant ids or chr_num and gene name and then the values are the ibd_program:filenames
    """
    # creating an empty dictionary to add values to 
    file_dict = dict()


    # iterate through the files to build the dictionary
    for ibd_file in ibd_file_list:

        print(ibd_file)

        chr_num: str = match_chr([r'.chr\d\d.', r'chr\d.'], 
                    ibd_file)

        # Finding the variant id of the file. file names are built so that the
        # variant id sits between the first "_" and the first "."
        shorten_ibd_file_string: str = get_ibd_file_substring(program_list, ibd_file)

        # TODO: This spot will break for the phenotype analysis
        identifier: str = determine_identifier(shorten_ibd_file_string, analysis_type)

        # using list comprehension to get all the files that contain that
        # variant and chromosome
        filter_ibd_file_list = [
            file for file in ibd_file_list
            if identifier in file and chr_num in file
        ]

        # match up the variants with the IBD program
        for ibd_program in program_list:

            # This goes through the three files in the filter_ibd_file_list
            for file in filter_ibd_file_list:

                # This checks to see if the variant id is in the dictionary
                # and that the ibd program is in the file
                if ibd_program in file and (
                        chr_num, identifier) not in file_dict.keys():

                    # If it is not then the
                    file_dict[(chr_num, identifier)] = set()

                    file_dict[(chr_num, identifier)].add("".join(
                        [ibd_program, ":", file]))

                elif ibd_program in file and (chr_num,
                                              identifier) in file_dict.keys():

                    file_dict[(chr_num, identifier)].add("".join(
                        [ibd_program, ":", file]))

    return file_dict