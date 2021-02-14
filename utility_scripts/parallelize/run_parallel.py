import multiprocessing as mp
from functools import partial
from dataclasses import dataclass
import utility_scripts


def parallelize_decorator(func):
    """decorator function that will be used to wrap the two methods for parallelization in the classes below.
    Parameter
    _________
    func : object
        function that the decorator will wrap around. This function should be a class method used 
        to parallelize the computation
    """

    def inner_func(self, *args):

        file_name: str = args[0]
        header_str: str = args[2]

        manager = mp.Manager()

        que = manager.Queue()
        pool = mp.Pool(int(self.workers))

        watcher = pool.apply_async(
            utility_scripts.listener,
            (que, "".join([self.output, file_name]), header_str))

        func(self,
             *args,
             que_object=que,
             pool_object=pool,
             manager_object=manager)

        que.put("kill")

        pool.close()

        pool.join()

    return inner_func


@dataclass
class Parallel_Runner:
    """parent class for code used to initial parallel runs

    Parameters
    __________
    workers : int
        number of cpu cores to be used during the computation

    output : str
        the path to the directory to output files into 
    """

    workers: int
    output: str


@dataclass
class Segment_Parallel_Runner(Parallel_Runner):
    """dataclass for the parallelization of the Shared_Segment_Formatter_CLI.py script

    Parameters
    __________
    ibd_format : str
        string listing the ibd program that was used. This will be either 
        hapibd or ilash

    min_CM : str
        string that lists the minimum centimorgan threshold to be used 
        throughout the computation

    file_list_dict : dict    
        dictionary containing a list of the iid_files, and a list of the 
        variant_info files. The keys in the dictionary will be "iid_files" 
        and "var_info_files"

    segment_file : str
        string that list the path to the segment file. This file should be 
        the output from either hapibd or ilash and should be a .match.gz or 
        a .ibd.gz
    """

    ibd_format: str
    min_CM: str
    file_list_dict: dict
    segment_file: str

    @parallelize_decorator
    def run_segments_parallel(self,
                              *args,
                              que_object=None,
                              pool_object=None,
                              manager_object=None):
        """function to run the computation in parallel
        Parameters
        __________
        file_name : str 
            string containing the filename for the output file

        parallel_func : object
            function that will be parallelized during the computation

        header_str : str
            string that will be the first row of the file at the file_name

        que_object : object
            que created by manager.Queue that waits for a string to be passed 
            to it and then passes that string to the listener function

        pool_object : object
            created by mp.Pool

        manager_object : object
            que manager created by mp.Manager

        """
        # expanding the second argument of the arg list into the parallel_func
        parallel_func: object = args[1]

        # expanding the file_list_dict variable into two separate list
        iid_file_list: list = self.file_list_dict["iid_files"]
        variant_info_list: list = self.file_list_dict["var_info_files"]

        func = partial(parallel_func, self.segment_file, self.output,
                       self.ibd_format, self.min_CM, iid_file_list, que_object)

        pool_object.map(func, variant_info_list)


@dataclass
class Haplotype_Parallel_Runner(Parallel_Runner):
    """child dataclass for the parallelization of the haplotype.py script

    Parameters
    __________
    file_list_dict : dict 
        dictionary that contains list of all the allpair files, the map files, 
        the carrier_files (The reformated single_variant_list.csv files), the
        ilash files, and the hapibd files.  This dictionary will have the keys 
        "allpair_files", "map_files", "carrier_files", "ilash_files", 
        "hapibd_files".

    network_file_path : str
        string that list the path to the network_groups.csv file

    variant_list : list 
        list of all the unique variants within the confirmed_carriers.txt files

    confirmed_carrier_file : str
        filepath to the confirmed_carriers.txt file

    """
    file_list_dict: dict
    network_file_path: str
    variant_list: list
    confirmed_carrier_file: str

    @parallelize_decorator
    def run_haplotypes_parallel(self,
                                *args,
                                que_object=None,
                                pool_object=None,
                                manager_object=None):
        """function to run the process in parallel
        Parameters
        __________
        file_name : str 
            string containing the filename for the output file

        parallel_func : object
            function that will be parallelized during the computation

        header_str : str
            string that will be the first row of the file at the file_name

        que_object : object
            que created by manager.Queue that waits for a string to be passed 
            to it and then passes that string to the listener function

        pool_object : object
            created by mp.Pool

        manager_object : object
            que manager created by mp.Manager

        """

        # starting the second que object
        variant_que = manager_object.Queue()

        variant_header: str = "variant\tchr\n"

        var_watcher = pool_object.apply_async(
            utility_scripts.listener, (variant_que, "".join([
                self.output, "nopairs_haplotype_analysis.txt"
            ]), variant_header))

        # expanding the file_list_dict to get all the
        allpair_file_list: list = self.file_list_dict["allpair_files"]
        map_file_list: list = self.file_list_dict["map_files"]
        carrier_file_list: list = self.file_list_dict["carrier_files"]
        ilash_file_list: list = self.file_list_dict["ilash_files"]
        hapibd_file_list: list = self.file_list_dict["hapibd_files"]

        # expanding the args[1] into the parallelized func
        parallel_func: object = args[1]

        func = partial(parallel_func, allpair_file_list, carrier_file_list,
                       map_file_list, ilash_file_list, hapibd_file_list,
                       que_object, variant_que, self.network_file_path,
                       self.confirmed_carrier_file)

        pool_object.map(func, self.variant_list)

        variant_que.put("kill")
