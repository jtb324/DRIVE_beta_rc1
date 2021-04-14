import pandas as pd
import os


def check_file_extension(file: str, file_extension_list: list) -> int:
    '''This function will see if the extension of the file is what you would expect.
    The function will return a 0 if it is not the expected file type or a one if it is'''

    input_file_extension: str = file.rfind(".")

    extension_handler_dict: dict = {
        True: 1,
        False: 0
    }

    return extension_handler_dict[input_file_extension in file_extension_list]
