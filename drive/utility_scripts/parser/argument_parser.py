# function that uses a generator to iterate through arguments passed to a function
def parser_generator(argument_dict: dict) -> tuple:
    """function that iterates through the key, value pairs in dictionary and 
        returns a tuple generator
    Parameters
    __________
    argument_dict: dict
        dictionary containing the key,value pairs of input arguments

    Returns
    _______
    tuple
        returns a tuple generator containing the key, value pair
    """
    for key, value in argument_dict.items():
        yield (key, value)
