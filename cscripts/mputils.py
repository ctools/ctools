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


# ==================== #
# Number of processes  #
# ==================== #
def nthreads(cls):
    """
    Determines the number of parallel processes to use.

    The number is based on the user parameter "nthreads" and the availability
    of the multiprocessing module.

    Parameters
    ----------
    cls : `~ctools.cscript`
        cscript class

    Returns
    -------
    nthreads : int
        Number of parallel processes
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