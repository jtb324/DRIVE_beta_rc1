from functools import wraps
import os


def check_dir_decorator(directory_name: str):
    def check_dir(func):
        """function to check if the directory exist"""
        @wraps(func)
        def inner_func(*args, **kwargs):
            
            output_path: str = args[0]["parameter_dict"]["output"]
            total_directory: str = os.path.join(output_path, directory_name)
            try:
                os.mkdir(total_directory)

            except FileExistsError:
                pass
            func(parameter_dict=args[0]["parameter_dict"])
        return inner_func
    return check_dir

def check_file_decorator(file_name: str):
    def check_file(func):
        """function to check if a file exist within a function"""
        @wraps(func)
        def inner_func(*args, **kwargs):

            output_path: str = kwargs["output"]
            total_filepath: str = os.path.join(output_path, file_name)
            # checking if the file path exist
            if os.path.exists(total_filepath):
                os.remove(total_filepath)

            func(*args, output=kwargs["output"])
        return inner_func
    return check_file
