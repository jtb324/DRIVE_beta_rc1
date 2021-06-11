def form_chr_num(number: int) -> str:
    """function to form the chromosome number for the key value in the get_file_dict dictionary
    Parameters
    __________
    number : int
        an integer number from the range function 
        
    Returns
    _______
    str
        returns a string of the form chrXX where X is a number
    """
    if number < 10:
        return "chr0"+str(number)
    else:
        return "chr"+str(number)

def find_match(file_list: list, chr_num: str, file_type: str) -> str:
    """Function to find the file that matches the chromosome number
    Parameters
    __________
    file_list : list 
        list of files from the get_file_dict function

    chr_num : str
        chromosome number that will be used to find the match in the 
        above files   

    file_type : str
        string that indicates if it is the map file, the carrier file, or the iid file
    Returns
    _______
    str
        returns the filepath of the matched file
    """

    # if the chromosome number is less than 10 you have to remove the
    # zero for the ibd files
    if int(chr_num[-2:]) < 10:
        alt_chr_num: str = chr_num.replace("0", "")
    else: 
        alt_chr_num = chr_num
    
    # using a dictionary to do the appropriate conditioning during the 
    # list comprehension for each file type
    file_handler = {
            "map": [file for file in file_list if "".join([chr_num, "_"]) in file],
            "ibd": [file for file in file_list if "".join([alt_chr_num, "."]) in file]
            }
    
    try:
        
        file: str = file_handler[file_type][0]

    except IndexError:

        return "None"
    
    return file

def get_file_dict(map_file_list: list,  ibd_file_list: list) -> dict:
    """Function to match up all of the files for the correct chromosome into a dictionary
    Parameters
    __________
    map_file_list : list
        a list of all the map files per chromosome as a result of running Plink
        
    ibd_file_list : list
        a list of all the files output by the specified ibd programs
        
    Returns
    _______
    dictionary
        returns a dictionary where the keys are the chromsomes and the values are dictionaries with the corresponding files
    """

    # forming a dictionary of dictionaries where the key is the 
    # chromosome number
    file_dict: dict = {form_chr_num(i):{} for i in range(1,22)}

    # get the files from each list
    for key in file_dict:

        file_dict[key].setdefault("map", find_match(map_file_list, key, "map"))
        file_dict[key].setdefault("ibd", find_match(ibd_file_list, key, "ibd"))
        
    return file_dict

def filter_empty_dictionaries(file_dict: dict) -> dict:
    """Function to filter the file dictionary to only those keys that contain a value
    Parameters
    __________
    file_dict : dict
        dictionary contain keys for each chromosome and dictionaries 
        with filepaths as values
    
    Returns
    _______
    dict
        returns the same dictionary as inputted but only the keys with 
        values not empty dictionaries
    """
    return {key:file_dict[key] for key in file_dict if "None" not in file_dict[key].values()}