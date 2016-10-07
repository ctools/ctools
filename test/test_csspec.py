#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csspec script.
#
# Copyright (C) 2016 Juergen Knoedlseder
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


# ============================ #
# Test class for csspec script #
# ============================ #
class Test(test):
    """
    Test class for csspec script

    This test class makes unit tests for the csspec script by using it
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
        self._inobs     = self._datadir + '/crab_cntmap_small.fits'
        self._expcube   = self._datadir + '/crab_expcube.fits'
        self._psfcube   = self._datadir + '/crab_psfcube.fits'
        self._edispcube = self._datadir + '/crab_edispcube.fits'
        self._bkgcube   = self._datadir + '/crab_bkgcube.fits'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csspec')

        # Append tests
        self.append(self._test_cmd, 'Test csspec on command line')
        self.append(self._test_python, 'Test csspec from Python')

        # Return
        return

    # Test csspec on command line
    def _test_cmd(self):
        """
        Test csspec on the command line
        """
        # Set script name
        csspec = self._script('csspec')

        # Setup csspec command
        cmd = csspec+' inobs="'+self._events+'"'+ \
                     ' inmodel="'+self._model+'"'+ \
                     ' srcname="Crab"'+ \
                     ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                     ' ebinalg="LOG" enumbins=5 emin=0.1 emax=100.0'+ \
                     ' outfile="csspec_cmd1.fits"'+ \
                     ' logfile="csspec_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('csspec_cmd1.fits', 5)

        # Setup csspec command
        cmd = csspec+' inobs="input_file_that_does_not_exist.fits"'+ \
                     ' inmodel="'+self._model+'"'+ \
                     ' srcname="Crab"'+ \
                     ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                     ' ebinalg="LOG" enumbins=5 emin=0.1 emax=100.0'+ \
                     ' outfile="csspec_cmd2.fits"'+ \
                     ' logfile="csspec_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csspec from Python
    def _test_python(self):
        """
        Test csspec from Python
        """
        # Set-up unbinned csspec
        spec = cscripts.csspec()
        spec['inobs']    = self._events
        spec['inmodel']  = self._model
        spec['srcname']  = 'Crab'
        spec['caldb']    = self._caldb
        spec['irf']      = self._irf
        spec['outfile']  = 'csspec_py1.fits'
        spec['ebinalg']  = 'LOG'
        spec['enumbins'] = 5
        spec['emin']     = 0.1
        spec['emax']     = 100.0
        spec['logfile']  = 'csspec_py1.log'
        spec['chatter']  = 2
        spec['publish']  = True

        # Run csspec script
        spec.logFileOpen()   # Make sure we get a log file
        spec.run()
        spec.save()

        # Check pull distribution file
        self._check_result_file('csspec_py1.fits', 5)

        # Set-up binned csspec
        spec = cscripts.csspec()
        spec['inobs']     = self._inobs
        spec['expcube']   = self._expcube
        spec['psfcube']   = self._psfcube
        spec['edispcube'] = self._edispcube
        spec['bkgcube']   = self._bkgcube
        spec['inmodel']   = self._model
        spec['srcname']   = 'Crab'
        spec['outfile']   = 'csspec_py2.fits'
        spec['ebinalg']   = 'LOG'
        spec['enumbins']  = 2
        spec['emin']      = 0.1
        spec['emax']      = 100.0
        spec['logfile']   = 'csspec_py2.log'
        spec['chatter']   = 3
        spec['publish']   = False

        # Execute csspec script
        spec.execute()

        # Check pull distribution file
        self._check_result_file('csspec_py2.fits', 2)

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, bins):
        """
        Check result file
        """
        # Open result file
        fits = gammalib.GFits(filename)

        # Get spectrum table
        spectrum = fits['SPECTRUM']

        # Check dimensions
        self.test_value(spectrum.nrows(), bins,
             'Check for %d rows in spectrum' % bins)
        self.test_value(spectrum.ncols(), 8,
             'Check for 8 columns in spectrum')

        # Return
        return
