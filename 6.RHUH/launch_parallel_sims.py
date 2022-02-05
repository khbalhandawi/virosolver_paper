import argparse
import numpy as np
import os
from multiprocessing import Pool
from itertools import repeat
import subprocess
from subprocess import PIPE,STDOUT

from typing import Any, Dict, Union, List, Dict, Callable

#=========================================================#
#                 Create empty directory                  #
#=========================================================#
def check_folder(folder,clear=False):
    """
    check if folder exists, make if not present

    Parameters
    ----------
    folder : str
        name of directory to check
    clear : bool, optional
        whether to clear the directory

    Raises
    ------
    Exception
        if folder could not be created or found
    """

    if not os.path.exists(folder):
        os.makedirs(folder)
    elif not os.path.exists(folder) and not clear:
        raise Exception('%s directory not found' %(folder))

    for f in os.listdir(folder):
        dirname = os.path.join(folder, f)
        if not os.path.isdir(dirname) and clear:
            os.remove(dirname)

#=========================================================#
#                  PARALLEL PROCESSING                    #
#=========================================================#
def starmap_with_kwargs(pool, fn, args_iter, kwargs_iter):
    """
    https://stackoverflow.com/questions/45718523/pass-kwargs-to-starmap-while-using-pool-in-python
    """
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(apply_args_and_kwargs, args_for_starmap)

def apply_args_and_kwargs(fn, args, kwargs):
    return fn(*args, **kwargs)

def parallel_sampling(func:Callable,
    vargs_iterator:List[List[Any]]=[]   ,vkwargs_iterator:List[Dict[str, Any]]=[],
    fargs:List[Any]=[]                  ,fkwargs:Dict[str, Any]={},
    num_threads:int=1):
    
    args_iter = []
    for vargs in vargs_iterator:
        args_iter += [vargs,*fargs]

    kwargs_iter = []
    for vkwargs in vkwargs_iterator:
        fkwargs.update(vkwargs)
        kwargs_iter += [fkwargs.copy()]

    if num_threads > 1:
        with Pool(num_threads) as pool:
            results = starmap_with_kwargs(pool, func, args_iter, kwargs_iter)
    else:
        results = []
        for args,kwargs in zip(args_iter,kwargs_iter):  
            result = func(*args, **kwargs)
            results.append(result)

    return results

#=========================================================#
#                      MCMC SCRIPT                        #
#=========================================================#
def system_command(command):
    """
    launches a string command in the windows terminal

    Parameters
    ----------
    command : str
        command to execute
    """
    #CREATE_NO_WINDOW = 0x08000000 # Create no console window flag
    p = subprocess.Popen(command,shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT,
                            ) # disable windows errors

    for line in iter(p.stdout.readline, b''):
        line = line.decode('utf-8')
        print(line.rstrip()) # print line by line
        # rstrip() to reomove \n separator

    streamdata = p.communicate()[0]
    rc = p.returncode
    # print('returncode : %i' %rc)
    if rc != 0:
        raise Exception('The following command failed : %s' %(command))

def run_mcmc(arguments,path,script,skip=False):
    """
    Executes an R script file with provided arguments

    Parameters
    ----------
    arguments : list
        list of argument and keyword pair
    path : str
        path to the main repository directory
    script : str
        path to the script file to be executed
    skip : bool, optional
        if True skips execution of the Abaqus script file, by default False
    """
    arguments_str = ' '.join(map(str,arguments)) # print variables as space delimited string
    script_file = os.path.join(path,script)

    results_folder = os.path.join(path,'..','mcmc_chains','5.real_leb_ct','sim_LEB_gp')
    check_folder(results_folder)
    out_log = os.path.join(results_folder,'a%i.Rout'%arguments[0])
    err_log = os.path.join(results_folder,'a%i.Rerr'%arguments[0])

    command = 'Rscript %s %s' %(script_file,arguments_str) + ' ' + '>' + out_log + ' ' + '2>' + err_log
    
    if not skip:
        system_command(command)

if __name__ == "__main__":

    N_CHAINS = 48

    vkwargs = []; i = 0
    for i in range(N_CHAINS):
        args_dict = {}
        args_dict['arguments']  = [i+1,]
        args_dict['path']       = os.getcwd()
        args_dict['script']     = '6.full_sim_recovery_cluster.R'
        vkwargs += [args_dict.copy()]

    vargs = []
    for i in range(len(vkwargs)):
        vargs += [[],]

    # set num_threads to one for serial processing
    results = parallel_sampling(run_mcmc,vargs_iterator=vargs,vkwargs_iterator=vkwargs,num_threads=8)
