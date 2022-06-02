# ==========================================================================
# Utility functions for multiprocessing handling
#
# Copyright (C) 2018 Luigi Tibaldo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import gammalib


# ===================== #
# Get number of threads #
# ===================== #
def nthreads(cls):
    """
    Determines the number of parallel processes to use

    The number is based on the user parameter "nthreads" and the availability
    of the multiprocessing module.

    Parameters
    ----------
    cls : `~ctools.cscript`
        cscript class

    Returns
    -------
    nthreads : int
        Number of threads
    """
    # Log multiprocessing configuration
    cls._log_header1(gammalib.NORMAL, 'Multiprocessing')

    # Try getting multiprocessing support
    try:
        from multiprocessing import cpu_count
        ncpus = cpu_count()
        cls._log_value(gammalib.NORMAL, 'Multiprocessing', 'available')
        cls._log_value(gammalib.EXPLICIT, 'CPUs available', ncpus)

        # Set processes to number of CPUs if requested by the user
        if cls['nthreads'].integer() == 0:
            nthreads = ncpus

        # ... otherwise use the number requested
        else:
            cls._log_value(gammalib.EXPLICIT, 'Processes requested',
                           cls['nthreads'].integer())
            nthreads = cls['nthreads'].integer()

    except:
        nthreads = 1
        cls._log_value(gammalib.NORMAL, 'Multiprocessing', 'not available')

    # Log number of processes
    cls._log_value(gammalib.NORMAL, 'Processes available', nthreads)

    # Return number of processes to be used
    return nthreads


# ============================ #
# Execute function in parallel #
# ============================ #
def process(nthreads, function, args):
    """
    Execute function in parallel

    Parameters
    ----------
    nthreads : int
        Number of parallel threads
    function : Python function
        Function to be executed in parallel
    args : list of tuples
        Function arguments for each evaluation

    Returns
    -------
    results : list of dicts
        List of function evaluation results
    """
    # Import Pool module from multiprocessing module
    from multiprocessing import Pool

    # Create a pool of processes
    pool = Pool(processes=nthreads)

    # Run function in parallel
    results = pool.map(function, args)

    # Close pool and join
    pool.close()
    pool.join()

    # Add CPU times of child processes to CPU time of parent
    nchilds = len(args)
    for i in range(nchilds):
        celapse = results[i][1]['celapse']
        args[i][0]._add_celapse(celapse)

    # Return results
    return results


# ======== #
# Function #
# ======== #
def mpfunc(args):
    """
    Multiprocessing function

    Parameters
    ----------
    args : tuple
        Tuple comprised of class, function name and function argument

    Returns
    -------
    result, info : tuple
        Tuple of function result and information
    """
    # Extract
    cls = args[0] # Class
    fct = args[1] # Class method name
    i   = args[2] # Function argument

    # Initialise thread logger
    cls._log.clear()
    cls._log.buffer_size(100000)

    # Compute function for argument
    cstart  = cls.celapse()
    result  = getattr(cls, fct)(i)
    celapse = cls.celapse() - cstart
    log     = cls._log.buffer()

    # Close logger
    cls._log.close()

    # Collect thread information
    info = {'celapse': celapse, 'log': log}

    # Return light curve bin result and thread information
    return result, info
