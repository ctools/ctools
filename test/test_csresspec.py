#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csresspec script.
#
# Copyright (C) 2017-2018 Luigi Tibaldo
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
import cscripts
from testing import test


# =============================== #
# Test class for csresspec script #
# =============================== #
class Test(test):
    """
    Test class for csresspec script

    This test class makes unit tests for the csresspec script by using it
    from the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set members
        self._inonoff     = self._datadir + '/onoff_obs.xml'
        self._onoff_model = self._datadir + '/onoff_model.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csresspec')

        # Append tests
        self.append(self._test_cmd, 'Test csresspec on command line')
        self.append(self._test_python, 'Test csresspec from Python')

        # Return
        return

    # Test csresspec on command line
    def _test_cmd(self):
        """
        Test csresspec on the command line.
        """
        # Set script name
        csresspec = self._script('csresspec')

        # Setup csresspec command
        cmd = csresspec + ' inobs="' + self._events + '"' + \
                          ' outfile="csresspec_cmd1.fits"' + \
                          ' inmodel="' + self._model + '"' + \
                          ' caldb="' + self._caldb + '" irf="' + self._irf + '"' + \
                          ' emin=1.0 emax=100.0 enumbins=2 ebinalg="LOG"' + \
                          ' mask=no algorithm="SUB"' + \
                          ' logfile="csresspec_cmd1.log" chatter=4'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
                         'Check successful execution from command line')

        # Check result file
        self._check_result_file('csresspec_cmd1.fits', 2, 5)

        # Setup csresspec command
        cmd = csresspec + ' inobs="event_file_that_does_not_exist.fits"' + \
                          ' outfile="csresspec_cmd2.fits"' + \
                          ' inmodel="' + self._model + '"' + \
                          ' caldb="' + self._caldb + '" irf="' + self._irf + '"' + \
                          ' emin=1.0 emax=100.0 enumbins=2 ebinalg="LOG"' + \
                          ' mask=no algorithm="SUBDIVSQRT"' + \
                          ' logfile="csresspec_cmd2.log" debug=yes chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check csresspec --help
        self._check_help(csresspec)

        # Return
        return

    # Test csresspec from Python
    def _test_python(self):
        """
        Test csresspec from Python
        """
        # Set-up csresspec for event list
        resspec = cscripts.csresspec()
        resspec['inobs']     = self._events
        resspec['outfile']   = 'csresspec_py1.fits'
        resspec['inmodel']   = self._model
        resspec['caldb']     = self._caldb
        resspec['irf']       = self._irf
        resspec['ebinalg']   = 'LOG'
        resspec['emin']      = 1.0
        resspec['emax']      = 100.0
        resspec['enumbins']  = 2
        resspec['mask']      = False
        resspec['algorithm'] = 'SUB'
        resspec['logfile']   = 'csresspec_py1.log'
        resspec['chatter']   = 2

        # Run csresspec script
        resspec.logFileOpen()  # Make sure we get a log file
        resspec.run()
        resspec.save()

        # Check result file
        self._check_result_file('csresspec_py1.fits', 2, 5)

        # Set-up csresspec for event list with computation of individual
        # model components
        resspec = cscripts.csresspec()
        resspec['inobs']      = self._events
        resspec['outfile']    = 'csresspec_py2.fits'
        resspec['inmodel']    = self._model
        resspec['caldb']      = self._caldb
        resspec['irf']        = self._irf
        resspec['ebinalg']   = 'LOG'
        resspec['emin']       = 1.0
        resspec['emax']       = 100.0
        resspec['enumbins']   = 2
        resspec['mask']       = False
        resspec['algorithm']  = 'SUBDIV'
        resspec['components'] = True
        resspec['logfile']    = 'csresspec_py2.log'
        resspec['chatter']    = 2

        # Run csresspec script
        resspec.execute()

        # Check result file
        self._check_result_file('csresspec_py2.fits', 2, 7)

        # Set-up csresspec for counts cube
        resspec = cscripts.csresspec()
        resspec['inobs']     = self._cntcube
        resspec['modcube']   = 'NONE'
        resspec['expcube']   = 'NONE'
        resspec['psfcube']   = 'NONE'
        resspec['edispcube'] = 'NONE'
        resspec['bkgcube']   = 'NONE'
        resspec['caldb']     = self._caldb
        resspec['irf']       = self._irf
        resspec['inmodel']   = self._model
        resspec['outfile']   = 'csresspec_py3.fits'
        resspec['mask']      = False
        resspec['algorithm'] = 'SUBDIVSQRT'
        resspec['logfile']   = 'csresspec_py3.log'
        resspec['chatter']   = 3

        # Run csresspec script
        resspec.execute()

        # Check result file
        self._check_result_file('csresspec_py3.fits', 20, 5)

        # Set-up csresspec for counts cube with source mask
        resspec = cscripts.csresspec()
        resspec['inobs']     = self._cntcube
        resspec['modcube']   = 'NONE'
        resspec['expcube']   = 'NONE'
        resspec['psfcube']   = 'NONE'
        resspec['edispcube'] = 'NONE'
        resspec['bkgcube']   = 'NONE'
        resspec['caldb']     = self._caldb
        resspec['irf']       = self._irf
        resspec['inmodel']   = self._model
        resspec['outfile']   = 'csresspec_py4.fits'
        resspec['mask']      = True
        resspec['ra']        = 83.63
        resspec['dec']       = 22.01
        resspec['rad']       = 0.2
        resspec['regfile']   = 'NONE'
        resspec['algorithm'] = 'SUBDIVSQRT'
        resspec['logfile']   = 'csresspec_py4.log'
        resspec['chatter']   = 3

        # Run csresspec script
        resspec.execute()

        # Check result file
        self._check_result_file('csresspec_py4.fits', 20, 5)

        # Set-up csresspec for On/Off observations
        resspec = cscripts.csresspec()
        resspec['inobs']     = self._inonoff
        resspec['inmodel']   = self._onoff_model
        resspec['outfile']   = 'csresspec_py5.fits'
        resspec['algorithm'] = 'SIGNIFICANCE'
        resspec['stack']     = False
        resspec['logfile']   = 'csresspec_py5.log'
        resspec['chatter']   = 3

        # Run csresspec script
        resspec.execute()

        # Check result file
        self._check_result_file('csresspec_py5.fits', 30, 8, nobs=2)

        # Set-up csresspec for On/Off observations
        # in stacked mode and with WSTAT background
        resspec = cscripts.csresspec()
        resspec['inobs']     = self._inonoff
        resspec['inmodel']   = self._onoff_model
        resspec['outfile']   = 'csresspec_py6.fits'
        resspec['algorithm'] = 'SIGNIFICANCE'
        resspec['statistic'] = 'WSTAT'
        resspec['stack']     = True
        resspec['logfile']   = 'csresspec_py6.log'
        resspec['chatter']   = 3

        # Run csresspec script
        resspec.execute()

        # Check result file
        self._check_result_file('csresspec_py6.fits', 30, 8)

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, enumbins, ncols, nobs=1):
        """
        Check result file
        """
        # Load residual file
        fits = gammalib.GFits(filename)

        # Check number of HDUs = observations (+ first empty image)
        self.test_value(fits.size(), 1 + nobs, 'Number of observations')

        # Check first hdu with residuals
        table = fits[1]

        # Check number of rows = number of energy bins
        self.test_value(table.nrows(), enumbins,
                        'Check for %d rows in residuals' % enumbins)

        # Check number of columns, 2 for energy bounds + 3 for 3D or 6 for
        # On/Off + 1 for each model component if components = yes
        self.test_value(table.ncols(), ncols,
                        'Check for %d columns in spectrum' % ncols)

        # Return
        return
