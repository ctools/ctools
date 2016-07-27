#! /usr/bin/env python
# ==========================================================================
# Performs unbinned maximum likelihood analysis for various models
#
# Copyright (C) 2011-2016 Juergen Knoedlseder
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
import sys
import ctools


# ================= #
# Analysis pipeline #
# ================= #
def pipeline(model_name):
    """
    Unbinned in-memory analysis pipeline

    This function implements an analysis pipeline that successively calls
    ctobssim, ctselect and ctlike without saving the intermediate results as
    FITS files on disk. All data is only hold in memory.
    """
    # Set script parameters
    caldb       = 'prod2'
    irf         = 'South_0.5h'
    ra          =   83.63
    dec         =   22.01
    rad_sim     =    5.0
    tstart      =    0.0
    tstop       =  180.0
    emin        =    0.1
    emax        =  100.0
    rad_select  =    3.0

    # Write model name
    print('*** Model: '+model_name+' ************************************')

    # Initialise timing
    wall_seconds = 0.0
    cpu_seconds  = 0.0

    # Simulate events
    sim = ctools.ctobssim()
    sim['inmodel'] = model_name
    sim['caldb']   = caldb
    sim['irf']     = irf
    sim['ra']      = ra
    sim['dec']     = dec
    sim['rad']     = rad_sim
    sim['tmin']    = tstart
    sim['tmax']    = tstop
    sim['emin']    = emin
    sim['emax']    = emax
    sim.run()

    # Select events
    select = ctools.ctselect(sim.obs())
    select['ra']   = ra
    select['dec']  = dec
    select['rad']  = rad_select
    select['tmin'] = tstart
    select['tmax'] = tstop
    select['emin'] = emin
    select['emax'] = emax
    select.run()

    # Perform maximum likelihood fitting
    like = ctools.ctlike(select.obs())
    like.run()

    # Show model fitting results
    print(like.obs().models()[0])

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':

    # Initialise flags
    need_help = False

    # Test for command line arguments
    print(sys.argv[0])
    if (len(sys.argv) > 1):
        if sys.argv[1] == '-h':
            need_help = True
        else:
            need_help = True

    # Print help if needed and exit
    if need_help:
        print('Usage: check_models.py [OPTIONS]')
        print('     -h       Display this usage message')
        sys.exit()

    # Dump header
    print('*******************************************')
    print('* Check models using an unbinned analysis *')
    print('*******************************************')

    # Perform analysis for Crab model
    pipeline('$CTOOLS/share/models/crab.xml')

    # Perform analysis for disk model
    pipeline('$CTOOLS/share/models/disk.xml')

    # Perform analysis for Gaussian model
    pipeline('$CTOOLS/share/models/gauss.xml')

    # Perform analysis for shell model
    pipeline('$CTOOLS/share/models/shell.xml')
