## ---------------------------
##
## Functions related to parallel processing
##
## Authors: Brian Gregor, Wenrui Li 
##          Boston University
##
## Date Created: 2020-07-31
##
## Email: bgregor@bu.edu
##
## ---------------------------
##
#  The sciris-based parallelization in covasim ran into memory errors
#  when run on 28-core machines so we wrote our own. This works without 
#  errors on 64-core machines and requires Python 3.8.x or newer due to
#  the use of the new multiprocessing.SharedMemory and 
#  multiprocessing.managers.SharedMemoryManager classes.  
#
#  By default this calls the simulation shrink() method which strips out the
#  People object to save RAM.  Reading in 1000 simulations with 30,000 people
#  for 109 days requires ~200 GB of RAM.
#
#  Read the docstring for the get_n_cores() function on using multiple cores.


import os
import multiprocessing as mp
from multiprocessing.managers import SharedMemoryManager

import tqdm
import sys
import pickle
import shutil
import string
import psutil 


# =============================================================================
#   Implement parallel simulation processing for covasim
#
#   Change 12/8/2020 - switched to a shared memory buffer for distributing
#    the original sim object to the Pool processes. 30% faster than the 
#    pickle file version previously used.
#
#   Change 7/22/2020 - write out simulation results as pickled files
#    during paralllel processing.  Then read the results into the master
#    process.  This saves gobs of RAM.  There's a moderate time penalty
#    due to dealing with the disk files.  Progress bars have been added.
# =============================================================================


__all__ = ['get_n_cores','parallel_run_sims']

#%%
def get_n_cores(use_physical_cores=False, cores_var='NSLOTS'):
    ''' Get the number of Cores that should be used.  This
    checks for an environment variable, NSLOTS, that is set
    using the SGE cluster queue software as used at BU's
    Shared Computing Cluster.  This should be modified to 
    fit your own needs, or you can set the NSLOTS variable
    manually before calling this code.  
    
    Optionally, set the use_physical_cores flag and this 
    will just use all available physical cores.
    '''
    # Note that cpu_count won't pick up the cores on a 2nd socket.  
    # If you know how many sockets/cores you have...just set the NSLOTS
    # variable.
    phys_cores = psutil.cpu_count(logical = False)
    if use_physical_cores:
        # The user knows best.
        return phys_cores
    ncores = 4 # polite limit on BU SCC login nodes.
    # but just to make sure this default isn't too many...
    if ncores > phys_cores:
        ncores = phys_cores
    # Environment checking for BU's SCC if this ran in the 
    # job queue.  
    if cores_var in os.environ:
        ncores = int(os.environ[cores_var])
    return ncores 
 


def parfun(sim, shrink=False):
    ''' Function used to evaluate simulations in parallel '''
    # Reset the random seed. MUST call this with 
    # a None to generate a new seed.
    sim.set_seed(None)
    sim.run()
    if shrink:
        sim.shrink() # remove the People to save RAM.
    return sim 


def unpickle_and_run_sim(args):
    ''' Run a simulation and pickle it to a file on the disk '''
    sim_shm, shrink = args
    # Reconstitute the sim
    sim = pickle.loads(sim_shm.buf)
    sim = parfun(sim, shrink)
    return sim


def gen_sims(sim_shm, n_runs, shrink):
    ''' Generator for the parallel sims'''
    for i in range(n_runs):
        yield sim_shm, shrink


def parallel_run_sims(sim,n_runs = 1, n_cores = 1, shrink = True, disable_progress = False):
    ''' Run n_runs simulations in parallel with the indicated 
        number of cores.  Returns a list of completed simulations.
        Use this in place of the MultiSim.run() method. 
        
        sim: a Covasim Sim object
        n_runs:  Number of simulations to run.  Default 1.
        n_cores: Number of cores to use. Default 1.
        shrink: Remove the People from completed simulations.  Saves RAM. Default True.
                This is always False if n_runs=1 for debugging.
        disable_progress: Disable the progress bars.  Default False.
        
        For n_runs=1 run it does not parallelize to help with debugging.'''
    sims_complete = []
    # Generator for the parallel sims.
    try:
        if n_runs > 1:
            print('Running %s simulations with %s cores.' % (n_runs, n_cores))
            sys.stdout.flush()
            # the context managers for the sharedmemorymanager and pool will
            # automatically clean up.
            with SharedMemoryManager() as smm:
                sim_str = pickle.dumps(sim)
                del sim  # save RAM
                sim_shm = smm.SharedMemory(size=len(sim_str))
                sim_shm.buf[:] = sim_str[:]
                del sim_str
                sims_data = gen_sims(sim_shm, n_runs, shrink)
                with mp.Pool(n_cores) as pool: # maxtasksperchild=1)
                    for res in tqdm.tqdm(pool.imap_unordered(unpickle_and_run_sim, sims_data), total=n_runs, disable=disable_progress):
                        sims_complete.append(res)
        else:
            # Don't run in parallel - this helps with debugging 1 process
            sims_complete.append(parfun(sim, shrink=False))
    except Exception as e:
        print('Failure running in parallel!')
        print(e)
        raise e
        
    print('Simulations complete.')
    return sims_complete
