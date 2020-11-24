import os
from file_creator_scripts.write_path import writePath


def check_dir(output_path: str, directory_name: str) -> str:
    network_directory = writePath(output_path, directory_name)

    try:
        os.mkdir(network_directory)
        print(f"Successfully created the {network_directory} directory")

    except FileExistsError:
        pass

    return network_directory
