from functools import wraps
from .readme_class import Readme

# This script will contain two decorators that is used to create READMEs for
# the given directories


def create_readmes(readme_text_list: list, output_path: str):
    """Function to create a readme
    Parameters
    __________
    readme_text_list : list
        list that contains all of the readme information

    output_path : str
        string that contains the output path to write the readme to
    """
    # creating a readme object using the Readme class
    readme = Readme("_README.md", output_path)
    # removing the previous file if it is there
    readme.rm_previous_file()
    # writing the output path to the log file
    readme.write_header(output_path)
    # creating a date tag
    readme.create_date_info()
    # adding all the text in the readme_text-list
    for readme_text in readme_text_list:
        readme.add_line(readme_text)


# decorator for class methods
def class_readme_generator(func):
    @wraps(func)
    def inner_func(self, **name):
        # getting the text from the
        readme_text: str = name["readme_txt"]
        output_path: str = name["output"]
        create_readmes(readme_text, output_path)
        func(self, **name)

    return inner_func


def func_readme_generator(func):
    @wraps(func)
    def inner_func(*args, **kwargs):
        
        readme_text: str = kwargs["parameter_dict"]["readme_text"]
        output_path: str = kwargs["parameter_dict"]["readme_output"]
        create_readmes(readme_text, output_path)

        if args:
            func(args[0])
        else:

            func(kwargs)

    return inner_func

