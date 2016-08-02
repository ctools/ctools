#! /usr/bin/env python
# ==========================================================================
# Generate TS distributions as function of energy for an OFF observation
#
# Copyright (C) 2011-2016 Jurgen Knodlseder
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
import cscripts


# ====================== #
# Create TS distribution #
# ====================== #
def create_ts(datadir, loge, emin, emax, ntrials=100, duration=180000.0,
              enumbins=0, log=False):
    """
    Create TS distribution.

    Parameters
    ----------
    datadir : str
        Data directory
    loge : float
        Logarithm of mean energy (TeV)
    emin : float
        Minimum energy (TeV)
    emax : float
        Maximum energy (TeV)
    ntrials : int, optional
        Number of MC samples
    durations : float, optional
        Observation duration (s)
    enumbins : int, optional
        Number of energy bins
    log : bool, optional
        Create log file(s)
    """
    # Generate output filename
    outfile = 'ts_'+str(loge)+'.dat'

    # Setup cstsdist tool
    tsdist = cscripts.cstsdist()
    tsdist['inmodel']  = datadir+'/crab.xml'
    tsdist['srcname']  = 'Crab'
    tsdist['outfile']  = outfile
    tsdist['ntrials']  = ntrials
    tsdist['caldb']    = 'prod2'
    tsdist['irf']      = 'South_50h'
    tsdist['ra']       = 83.63
    tsdist['dec']      = 22.01
    tsdist['emin']     = emin
    tsdist['emax']     = emax
    tsdist['enumbins'] = enumbins
    tsdist['tmin']     = 0
    tsdist['tmax']     = duration
    tsdist['rad']      = 5.0
    tsdist['npix']     = 200
    tsdist['binsz']    = 0.05

    # Optionally open the log file
    if log:
        tsdist.logFileOpen()

    # Run tool
    tsdist.run()

    # Return
    return


# ==================== #
# Make TS distribution #
# ==================== #
def make_ts_distribution():
    """
    Make TS distribution
    """
    # Set usage string
    usage = 'make_ts_distributions.py [-n ntrials] [-e enumbins]'\
            ' [-t duration] [-m max_threads] [-d datadir]'

    # Set default options
    options = [{'option': '-n', 'value': '100'},
               {'option': '-e', 'value': '0'},
               {'option': '-t', 'value': '180000.0'},
               {'option': '-m', 'value': '1'},
               {'option': '-d', 'value': 'data'}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    ntrials     = int(options[0]['value'])
    enumbins    = int(options[1]['value'])
    duration    = float(options[2]['value'])
    max_threads = int(options[3]['value'])
    datadir     = options[4]['value']

    # Loop over energy bands. The energy bands are those that are also
    # used for sensitivity computation.
    for ieng in range(20):

        # Set energies
        loge  = -1.7 + ieng * 0.2
        emean = pow(10.0, loge)
        emin  = pow(10.0, loge-0.1)
        emax  = pow(10.0, loge+0.1)
        if loge < 0:
            loge = 'm'+str(abs(loge))
        else:
            loge = 'p'+str(abs(loge))

        # Create TS
        create_ts(datadir, loge, emin, emax, ntrials=ntrials, enumbins=enumbins,
                  duration=duration)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Make TS distribution
    make_ts_distribution()
