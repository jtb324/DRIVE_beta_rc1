# This file will contain the parent class for generating readme and log files
import os
from os import path
import getpass
from datetime import datetime

# This is the base class for creating documentation


class Documentation:
    '''This class will be the parent class for the documentation throughout
    the program'''
    def __init__(self, file_name: str, output_path: str):

        self.file_name = "".join([output_path, file_name])
        current_time = datetime.now()
        self.day = current_time.strftime("%b %d, %Y")
        self.time = current_time.strftime("%H:%M:%S")
        self.user = getpass.getuser()

    def rm_previous_file(self):
        '''This function will remove the readme from a previous run'''

        if path.exists(self.file_name):

            os.remove(self.file_name)

    def create_date_info(self):
        '''This function will create a file with a header'''

        with open(self.file_name, "a+") as file:

            # Writing header information to the file
            file.write(
                "################################################################\n\n"
            )
            file.write(
                "DRIVE: Distant Relatedness for Identification and Variant Evaluation\n"
            )
            file.write(f"run started on {self.day} at {self.time}\n")
            file.write(f"run started by user: {self.user}\n\n")
            file.write(
                "################################################################\n\n"
            )


class Readme(Documentation):
    '''This class will extend the Documentation class'''
    def __init__(self, file_name: str, output_path: str):

        super().__init__(file_name, output_path)

    def write_header(self, directory_name: str):
        '''This function will write a fairly genetic header for the file'''

        with open(self.file_name, "a+") as readme_file:

            readme_file.write(
                f"# README for the {directory_name} directory:\n\n")

    def add_line(self, info_str: str):
        '''This function will the info_str to a new line in the document'''

        with open(
                self.file_name,
                "a+",
        ) as readme_file:
            readme_file.write(info_str + "\n")


class LogFile(Documentation):
    '''This class is used to create log files and is extended the
    parent documentation class'''
    def __init__(self, file_name: str, output_path: str):

        super().__init__(file_name, output_path)

    def write_header(self):
        '''This function creates the header line for the log file'''

        with open(self.file_name, "a+") as readme_file:

            readme_file.write("Log file for the DRIVE program:\n")

    @staticmethod
    def create_format_string(log_tag: str, text: str, time: str) -> str:
        '''This function is responsible for '''

        return f"{time} -- {log_tag}:   {text}\n"

    def write_to_file(self, log_str: str):
        '''This function will actually write the log string to a file'''

        with open(self.file_name, "a+") as log_file:

            log_file.write(log_str)

    def add_newline(self, log_tag: str, text: str):
        '''This will be the main function that is used to add a new line to the
        log file'''

        # Getting the time when the line will be added to the logfile
        datetime_obj = datetime.now()
        time = datetime_obj.strftime("%H:%M:%S")

        log_str: str = self.create_format_string(log_tag, text, time)

        self.write_to_file(log_str)
