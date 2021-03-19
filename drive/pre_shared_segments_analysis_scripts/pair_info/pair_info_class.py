from dataclasses import dataclass


class Pair_Info_Class:
    """class to organize information about the pair 1 and pair 2 
    Parameters
    __________
    pair_string : str
        string containing the grid ID for Pair 1 and pair 2

    variant_id : str
        the variant that the pairs are being tested for
    
    chromo_num : str
        string describing the chromosome number that the shared 
        segment is on
    allpair_filepath : str
        string containing the filepath that the allpair file will be output to
    """

    def __init__(self, pair_string: str, variant_id: str, chromo_num: str, allpair_filepath: str):
        self.pair_1, self.pair_2 = self.split_pair_str(pair_string)
        self.variant_id: str = variant_id
        self.chromo_num: str = chromo_num
        self.output_path: str = allpair_filepath
    
    @staticmethod
    def split_pair_str(pair_string) -> list:
        """Function to split the string of pairs 
        Parameters
        __________
        pair_string : str
            string containing the grid ids for pair 1 and pair 2
        Returns
        _______
        list 
            list containing the id for pair 1 and pair 2
        """
        # getting the string of just the two pairs
        pairs_str: str = pair_string.split("\t")[4]

        # splitting the above string into a list with pair 1 and
        # pair 2
        pairs_list: list = pairs_str.split(" ")

        return pairs_list