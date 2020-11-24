#######################################################
# Import modules
import csv

# Import functions from other files
import file_creator_scripts

#######################################################
# Creating the class to write a dictionary to csv file


class Csv_Writer_Object:

    def __init__(self, var_dictionary, write_directory, file_name, logger):
        '''This is invoked when the class is instantiated. The class requires a dictionary to be passed. A name for the directory it will be written to, a filename for the output file and then the logger so that the program can record where the file is being output to in a log file.'''

        self.var_dict = var_dictionary
        self.log_file = logger
        self.output_directory = write_directory
        self.output_file_name = file_name

    def log_file_path(self):
        '''This function writes a message to a log file indicating where the output file can be found'''
        self.log_file.info('Writing the dictionary to the file path {}'.format(
            file_creator_scripts.writePath(self.output_directory, self.output_file_name)))

    def write_to_csv(self):
        '''This function writes the passed dictionary and writes it to a csv file.'''

        with open(file_creator_scripts.writePath(self.output_directory, self.output_file_name), 'w') as csvfile:

            # this line creates a writer object that the multiVarDict will be written to
            writer = csv.writer(csvfile)

            for key, value in self.var_dict.items():  # This line iterates through the keys and values in the multiVarDict

                # This line writes the keys and values to a row
                writer.writerow([key, value])

        csvfile.close()  # This closes the csv file
