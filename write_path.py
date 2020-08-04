import os


def writePath(write_location, file_name):
    '''This section creates a function that creates a write path that can be used. This helps keep the code DRY'''

    # This line just joins the path of the directory to write to and the filename for a complete string.
    total_var_directory = os.path.join(
        write_location, file_name)

    return total_var_directory
