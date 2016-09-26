#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cttsmap tool
#
# Copyright (C) 2014-2016 Juergen Knoedlseder
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
import os
import gammalib
import ctools
from testing import test


# =========================== #
# Test class for cttsmap tool #
# =========================== #
class Test(test):
    """
    Test class for cttsmap tool
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('cttsmap')

        # Append tests
        self.append(self._test_cmd, 'Test cttsmap on command line')
        self.append(self._test_python, 'Test cttsmap from Python')

        # Return
        return

    # Test cttsmap on command line
    def _test_cmd(self):
        """
        Test cttsmap on the command line
        """
        # Set tool name
        cttsmap = self._tool('cttsmap')

        # Setup cttsmap command
        cmd = cttsmap+' inobs="'+self._events+'"'+ \
                      ' inmodel="'+self._model+'" srcname="Crab"'+ \
                      ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                      ' outmap="cttsmap_cmd1.fits"'+ \
                      ' nxpix=3 nypix=3 binsz=0.02'+ \
                      ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                      ' logfile="cttsmap_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('cttsmap_cmd1.fits')

        # Setup cttsmap command
        cmd = cttsmap+' inobs="event_file_that_does_not_exist.fits"'+ \
                      ' inmodel="'+self._model+'" srcname="Crab"'+ \
                      ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                      ' outmap="cttsmap_cmd2.fits"'+ \
                      ' nxpix=3 nypix=3 binsz=0.02'+ \
                      ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                      ' logfile="cttsmap_cmd2.log" chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup cttsmap command
        cmd = cttsmap+' inobs="'+self._events+'"'+ \
                      ' inmodel="'+self._model+'" srcname="Crab"'+ \
                      ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                      ' errors=yes'+ \
                      ' outmap="cttsmap_cmd3.fits"'+ \
                      ' nxpix=3 nypix=3 binsz=0.02'+ \
                      ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                      ' logfile="cttsmap_cmd3.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('cttsmap_cmd3.fits')

        # Setup cttsmap --help option
        cmd = cttsmap+' --help'

        # Check if execution was successful in case that the CTOOLS
        # environment variable was set or failed otherwise
        if 'CTOOLS' in os.environ:
            self.test_value(self._execute(cmd), 0,
                 'Check successful execution with --help option')
        else:
            self.test_assert(self._execute(cmd) != 0,
                 'Check execution failure with --help option')

        # Return
        return

    # Test cttsmap from Python
    def _test_python(self):
        """
        Test cttsmap from Python
        """
        # Allocate empty cttsmap tool
        tsmap = ctools.cttsmap()

        # Check that empty cttsmap tool holds an empty observation and an
        # empty TS map
        self._check_obs(tsmap.obs(), nobs=0)
        self._check_ts_map(tsmap.tsmap(), nx=0, ny=0)

        # Check that saving does not nothing
        tsmap['outmap']  = 'cttsmap_py0.fits'
        tsmap['logfile'] = 'cttsmap_py0.log'
        tsmap.logFileOpen()
        tsmap.save()
        self.test_assert(not os.path.isfile('cttsmap_py0.fits'),
             'Check that no TS map has been created')

        # Check that clearing does not lead to an exception or segfault
        tsmap.clear()

        # Now set cttsmap parameters
        tsmap['inobs']    = self._events
        tsmap['inmodel']  = self._model
        tsmap['srcname']  = 'Crab'
        tsmap['caldb']    = self._caldb
        tsmap['irf']      = self._irf
        tsmap['nxpix']    = 3
        tsmap['nypix']    = 3
        tsmap['binsz']    = 0.02
        tsmap['coordsys'] = 'CEL'
        tsmap['proj']     = 'CAR'
        tsmap['xref']     = 83.63
        tsmap['yref']     = 22.01
        tsmap['outmap']   = 'cttsmap_py1.fits'
        tsmap['logfile']  = 'cttsmap_py1.log'
        tsmap['chatter']  = 2

        # Run cttsmap tool
        tsmap.logFileOpen()   # Make sure we get a log file
        tsmap.run()
        tsmap.save()

        # Check result file
        self._check_result_file('cttsmap_py1.fits')

        # Copy cttsmap tool
        cpy_tsmap = tsmap.copy()

        # Check observation container and TS map of copy
        self._check_obs(cpy_tsmap.obs())
        self._check_ts_map(cpy_tsmap.tsmap())

        # Execute copy of cttsmap tool again, now with a higher chatter
        # level than before
        cpy_tsmap['binmin']  = 0
        cpy_tsmap['binmax']  = 1
        cpy_tsmap['outmap']  = 'cttsmap_py2.fits'
        cpy_tsmap['logfile'] = 'cttsmap_py2.log'
        cpy_tsmap['chatter'] = 3
        cpy_tsmap['publish'] = True
        cpy_tsmap.logFileOpen()  # Needed to get a new log file
        cpy_tsmap.execute()

        # Check result file
        self._check_result_file('cttsmap_py2.fits')

        # Now clear copy of cttsmap tool
        cpy_tsmap.clear()

        # Check that the cleared copy has also cleared the observation
        # container and TS map
        self._check_obs(cpy_tsmap.obs(), nobs=0)
        self._check_ts_map(cpy_tsmap.tsmap(), nx=0, ny=0)

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nx=3, ny=3):
        """
        Check result file

        Parameters
        ----------
        filename : str
            TS map file name
        nx : int, optional
            Number of X pixels
        ny : int, optional
            Number of Y pixels
        """
        # Open result file
        fits = gammalib.GFits(filename)

        # Get HDUs
        ts    = fits['Primary']
        pre   = fits['Prefactor']
        index = fits['Index']

        # Check dimensions
        self.test_value(ts.naxis(), 2, 'Check for TS map dimensions')
        self.test_value(ts.naxes(0), nx, 'Check for TS map pixels in X')
        self.test_value(ts.naxes(1), ny, 'Check for TS map pixels in Y')
        self.test_value(pre.naxis(), 2, 'Check for Prefactor map dimensions')
        self.test_value(pre.naxes(0), nx, 'Check for Prefactor map pixels in X')
        self.test_value(pre.naxes(1), ny, 'Check for Prefactor map pixels in Y')
        self.test_value(index.naxis(), 2, 'Check for Index map dimensions')
        self.test_value(index.naxes(0), nx, 'Check for Index map pixels in X')
        self.test_value(index.naxes(1), ny, 'Check for Index map pixels in Y')

        # Return
        return

    # Check observation
    def _check_obs(self, obs, nobs=1):
        """
        Check observation

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Models
        nobs : int, optional
            Expected number of observations
        """
        # Check number of observations
        self.test_value(obs.size(), nobs, 'Check number of observations')

        # Return
        return

    # Check TS map
    def _check_ts_map(self, tsmap, nx=3, ny=3):
        """
        Check TS map

        Parameters
        ----------
        tsmap : `~gammalib.GSkyMap`
            TS map
        nx : int, optional
            Number of X pixels
        ny : int, optional
            Number of Y pixels
        """
        # Determine number of maps
        if nx > 0 and ny > 0:
            nmaps = 1
        else:
            nmaps = 0

        # Check dimensions
        self.test_value(tsmap.nmaps(), nmaps, 'Check number of maps')
        self.test_value(tsmap.nx(), nx, 'Check for number of X pixels')
        self.test_value(tsmap.ny(), ny, 'Check for number of Y pixels')

        # Return
        return
