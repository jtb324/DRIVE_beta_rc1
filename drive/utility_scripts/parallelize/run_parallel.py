import multiprocessing as mp
from functools import partial
from dataclasses import dataclass
from typing import List, Dict
import utility_scripts


# def parallelize_decorator(func):
#     """decorator function that will be used to wrap the two methods for parallelization in the classes below.
#     Parameter
#     _________
#     func : object
#         function that the decorator will wrap around. This function should be a class method used 
#         to parallelize the computation
#     """
#     def inner_func(self, *args):

#         file_name: str = args[0]
#         header_str: str = args[2]

#         manager = mp.Manager()

#         que = manager.Queue()
#         pool = mp.Pool(int(self.workers))

#         watcher = pool.apply_async(
#             utility_scripts.listener,
#             (que, "".join([self.output, file_name]), header_str))

#         func(self,
#              *args,
#              que_object=que,
#              pool_object=pool,
#              manager_object=manager)

#         que.put("kill")

#         pool.close()

#         pool.join()

#     return inner_func


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
        string listing the ibd program that was used. This will be either hapibd or ilash

    min_CM : str
        string that lists the minimum centimorgan threshold to be used throughout the computation

    file_list_dict : dict    
        dictionary of dictionaries where the outer key is the 
        variant and the inner dictionary has two keys, base_pos and iid_list which have the variants base position and 
        the list of grids that carry the variant according to the
        mega array

    segment_file : str
        string that list the path to the segment file. This file should be the output from either hapibd or ilash and should be a .match.gz or a .ibd.gz
    """

    ibd_format: str
    min_CM: str
    file_list_dict: Dict[str, Dict]
    segment_file: str

    def run_segments_parallel(self,
        parallel_func: object, 
        file_name: str,
        header: str) -> Dict:
        """function to run the computation in parallel
        Parameters
        __________
        file_name : str 
            string containing the filename for the output file

        parallel_func : object
            function that will be parallelized during the computation

        header_str : str
            string that will be the first row of the file at the file_name

        """
        # expanding the second argument of the arg list into the parallel_func
        manager = mp.Manager()

        que = manager.Queue()

        pool = mp.Pool(int(self.workers))

        pair_info_dict: Dict[str, Dict] = manager.dict()
        
        # get all the variants from the file_list_dict
        variant_list: list = self.file_list_dict.keys()
        
        _ = pool.apply_async(
            utility_scripts.listener,
            (que, "".join([self.output, file_name]), header))

        # func = partial(parallel_func, self.segment_file, self.output,
        #                self.ibd_format, self.min_CM, self.file_list_dict, pair_info_dict, que)
        
        for variant in variant_list:

            pool.apply_async(parallel_func, args=(variant, self.segment_file, self.output,
                       self.ibd_format, self.min_CM, self.file_list_dict, pair_info_dict, que))
                       
        # pool.map(func, variant_list)

        que.put("kill")

        pool.close()

        pool.join()
        
        return pair_info_dict
        

def parallelize_test(func: object, cpu_count: int, pair_info_dict: Dict[str, Dict], combined_info_list: List,  output: str = None, que_object: bool = False):

        pool = mp.Pool(cpu_count)

        func_partial = partial(func, pair_info_dict)
        
        pool.map(func_partial, combined_info_list)
        
        pool.close()

        pool.join()
