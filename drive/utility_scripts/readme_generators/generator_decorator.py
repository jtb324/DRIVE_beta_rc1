from functools import wraps
from .readme_class import Readme

# This script will contain two decorators that is used to create READMEs for
# the given directories


def create_readmes(readme_text: str, output_path: str):
    readme = Readme("_README.md", output_path)
    readme.rm_previous_file()
    readme.write_header(output_path)
    readme.create_date_info()
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
        readme_text: str = kwargs["readme_text"]
        output_path: str = kwargs["readme_output"]
        create_readmes(readme_text, output_path)
        func(*args, **kwargs)

    return inner_func
