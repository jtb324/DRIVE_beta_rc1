from dataclasses import dataclass
import utility_scripts

@dataclass
class Input_Files:
    """Class for all of the input files 
    Parameters
    __________
    carrier_file_path : str
        string listing the filepath to the directory containing the .single_variant_carrier.csv files
    
    hapibd_file : str
        string listing the filepath to the hapibd file that will be used
    
    ilash_file : str
        string listing the filepath to the iLash file that will be used

    map_file : str
        string that will list the filepath to the map file that will be used
    """
    carrier_file_path : str
    hapibd_file : str
    ilash_file : str
    map_file : str

    def get_carrier_file_list(self) -> list:
        """Function to get a list of single_variant_carrier.csv
        Returns
        _______
        list
            returns a list of all the single_variant_carriers.csv file
        """
        car_file_list: list = utility_scripts.get_file_list(self.carrier_file_path, "*.csv")

        return car_file_list