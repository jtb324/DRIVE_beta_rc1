from functools import wraps
from .readme_class import Readme

# This script will contain two decorators that is used to create READMEs for
# the given directories


# decorator for class methods
def class_readme_generator(func):
    @wraps(func)
    def inner_func(self, **name):
        # getting the text from the
        readme_text: str = name["readme_txt"]
        readme = Readme("_README.md",
                        "".join([self.output, "plink_output_files/"]))
        readme.rm_previous_file()
        readme.write_header("".join(["plink_output_files/"]))
        readme.create_date_info()
        readme.add_line(readme_text)
        func(self, **name)

    return inner_func
