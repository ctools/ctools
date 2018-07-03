# ==========================================================================
# Utility functions for multiprocessing handling
#
# Copyright (C) 2018- Luigi Tibaldo
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
    Determines the number of parallel processes to use,
    based on user parameters and availability of multiprocessing module

    :param cls: `~ctools.cscript`
        cscript class
    :return: int
        number of parallel processes
    """

    # Log multiprocessing configuration
    cls._log_header1(gammalib.NORMAL, "Multiprocessing")

    try:
        from multiprocessing import cpu_count
        ncpus = cpu_count()
        cls._log_value(gammalib.EXPLICIT, 'Number of CPUs available', ncpus)

        # Set processes to number of CPUs if requested by the user
        if cls['nthreads'].integer() == 0:
            nthr = ncpus

        # Otherwise use the number requested
        else:
            cls._log_value(gammalib.EXPLICIT, 'Number of processes requested by the user',
                           cls['nthreads'].integer())
            nthr = cls['nthreads'].integer()

    except:
        nthr = 1
        cls._log_value(gammalib.EXPLICIT,
                       'Multiprocessing not available. Number of processes', nthr)


    # Log number of processes
    cls._log_value(gammalib.NORMAL, 'Number of parallel processes in use', nthr)

    # return number of processes to be used
    return nthr
