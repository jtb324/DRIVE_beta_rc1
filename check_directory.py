import os
from write_path import writePath


def check_dir(output_path, directory_name):
    network_directory = writePath(output_path, directory_name)

    try:
        os.mkdir(network_directory)
        print("Successfully created the {} directory".format(network_directory))

    except FileExistsError:
        print("{} is already an existing directory".format(network_directory))

    return network_directory
