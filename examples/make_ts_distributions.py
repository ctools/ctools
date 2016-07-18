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


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Get input arguments
    usage = 'make_ts_distributions.py [-n ntrials] [-e enumbins] [-m max_threads]'
    if len(sys.argv) < 1:
        print(usage)
        sys.exit()

    # Set default parameters
    ntrials     = 100
    enumbins    = 0
    duration    = 180000.0
    max_threads = 1
    datadir     = 'data'

    # Parameter dictionnary
    pars = [{'option': '-n', 'value': ntrials},
            {'option': '-e', 'value': enumbins},
            {'option': '-d', 'value': duration},
            {'option': '-m', 'value': max_threads},
            {'option': '-p', 'value': datadir}]

    # Gather parameters from command line
    i = 1
    while i < len(sys.argv):

        # Search for options
        for par in pars:
            if sys.argv[i] == par['option']:
                if len(sys.argv) > i+1:
                    i += 1
                    try:
                        par['value'] = int(sys.argv[i])
                    except:
                        print(usage)
                        sys.exit()
                else:
                    print(usage)
                    sys.exit()

        # Next item
        i += 1

    # Recover parameters
    ntrials     = pars[0]['value']
    enumbins    = pars[1]['value']
    duration    = pars[2]['value']
    max_threads = pars[3]['value']
    datadir     = pars[4]['value']

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
