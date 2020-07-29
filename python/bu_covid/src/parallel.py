# Functions related to parallel processing

import os
import multiprocessing as mp
import tempfile
import tqdm
import sys
import pickle
import shutil
import secrets
import string
import psutil 


# =============================================================================
#   Implement parallel simulation processing for covasim
#   Change 7/22/2020 - write out simulation results as pickled files
#    during paralllel processing.  Then read the results into the master
#    process.  This saves gobs of RAM.  There's a moderate time penalty
#    due to dealing with the disk files.  Progress bars have been added.
# =============================================================================


__all__ = ['get_n_cores','parallel_run_sims']

#%%
def get_n_cores(use_physical_cores=False):
    ''' Get the number of Cores that should be used.  This
    checks for an environment variable, NSLOTS, that is set
    using the SGE cluster queue software as used at BU's
    Shared Computing Cluster.  This should be modified to 
    fit your own needs, or you can set the NSLOTS variable
    manually before calling this code.  
    
    Optionally, set the use_physical_cores flag and this 
    will just use all available physical cores.
    '''
    # Note that this won't pick up the cores on a 2nd socket.  If you know 
    # how many sockets/cores you have...just set the NSLOTS
    # variable.
    phys_cores = psutil.cpu_count(logical = False)
    if use_physical_cores:
        # The user knows best.
        return phys_cores
    ncores = 4 # polite limit on BU SCC login nodes.
    # but just to make sure this default isn't too many...
    if ncores > phys_cores:
        ncores = phys_cores
    # Environment checking for BU's SCC
    if 'NSLOTS' in os.environ:
        ncores = int(os.environ['NSLOTS'])
    return ncores 



def get_random_string(length):
    ''' Get a random string '''
    secure_str = ''.join((secrets.choice(string.ascii_letters) for i in range(length)))
    return secure_str


def parfun(sim):
    ''' Function used to evaluate simulations in parallel '''
    # Reset the random seed. MUST call this with 
    # a None to generate a new seed.
    sim.set_seed(None)
    sim.run()
    return sim 


def run_and_pickle_sim(args):
    ''' Run a simulation and pickle it to a file on the disk '''
    sim, tmpdir, n, shrink = args
    sim = parfun(sim)
    if shrink:
        sim.shrink() # remove the People to save RAM.
    out_file = os.path.join(tmpdir,'sim_%s.pkl' % n)
    with open(out_file,'wb') as f:
        pickle.dump(sim, f, protocol = pickle.HIGHEST_PROTOCOL)
    return out_file


def unpickle_sims(sims_paths, disable_progress = False):
    ''' Unpickle the sims and delete the pickled files '''
    sims_complete = []
    print('Unpickling completed simulations.')
    sys.stdout.flush()
    for sim_pkl in tqdm.tqdm(sims_paths, disable = disable_progress):
        try:
            with open(sim_pkl,'rb') as f:
                sim = pickle.load(f)
                sims_complete.append(sim)
            os.unlink(sim_pkl)
        except Exception as e:
            print('Error unpickling %s' % sim_pkl)
            print(e)
    return sims_complete


def gen_sims(sim, tmpdir, n_runs, shrink):
    ''' Generator for the parallel sims'''
    for i in range(n_runs):
        yield sim, tmpdir, i, shrink


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
    tmpdir = os.path.join(tempfile.gettempdir(), get_random_string(16))
    os.makedirs(tmpdir, exist_ok = True)
    sims_complete = []
    # Generator for the parallel sims.
    try:    
        sims_data = gen_sims(sim, tmpdir, n_runs, shrink)
        if n_runs > 1:
            print('Running %s simulations with %s cores.' % (n_runs, n_cores))
            sys.stdout.flush()
            # Start up the pool.  Each Python proc runs exactly 1 simulation and then
            # gets shut down to guarantee there's no weirdness with shared lists.
            pool = mp.Pool(n_cores,maxtasksperchild=1)
            sims_paths = []
            for res in tqdm.tqdm(pool.imap_unordered(run_and_pickle_sim, sims_data), total=n_runs, disable=disable_progress):
                sims_paths.append(res)
            pool.close()
            pool.join()
        else:
            # Don't run in parallel - this helps with debugging 1 process
            sims_paths = [run_and_pickle_sim( (sim, tmpdir, 0, False) )]  
        # Gather up the completed sims
        sims_complete = unpickle_sims(sims_paths, disable_progress)
    except Exception as e:
        print('Failure running in parallel!')
        print(e)
        raise e
    finally:
        # Guaranteed cleanup. Delete tmpdir.
        shutil.rmtree(tmpdir)  
        
    print('Simulations complete.')
    return sims_complete
