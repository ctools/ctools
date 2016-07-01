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
        cmd = csspec+' inobs="data/crab_events.fits"'+ \
                     ' inmodel="data/crab.xml"'+ \
                     ' srcname="Crab" caldb="prod2" irf="South_0.5h"' + \
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
                     ' inmodel="data/crab.xml"'+ \
                     ' srcname="Crab" caldb="prod2" irf="South_0.5h"' + \
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
        spec['inobs']    = 'data/crab_events.fits'
        spec['inmodel']  = 'data/crab.xml'
        spec['srcname']  = 'Crab'
        spec['caldb']    = 'prod2'
        spec['irf']      = 'South_0.5h'
        spec['outfile']  = 'csspec_py1.fits'
        spec['ebinalg']  = 'LOG'
        spec['enumbins'] = 5
        spec['emin']     = 0.1
        spec['emax']     = 100.0
        spec['logfile']  = 'csspec_py1.log'
        spec['chatter']  = 2

        # Run csspec script
        spec.logFileOpen()   # Make sure we get a log file
        spec.run()
        spec.save()

        # Check pull distribution file
        self._check_result_file('csspec_py1.fits', 5)

        # Set-up binned csspec
        spec = cscripts.csspec()
        spec['inobs']     = 'data/crab_cntmap_small.fits'
        spec['expcube']   = 'data/crab_expcube.fits'
        spec['psfcube']   = 'data/crab_psfcube.fits'
        spec['edispcube'] = 'data/crab_edispcube.fits'
        spec['bkgcube']   = 'data/crab_bkgcube.fits'
        spec['inmodel']   = 'data/crab.xml'
        spec['srcname']   = 'Crab'
        spec['outfile']   = 'csspec_py2.fits'
        spec['ebinalg']   = 'LOG'
        spec['enumbins']  = 2
        spec['emin']      = 0.1
        spec['emax']      = 100.0
        spec['logfile']   = 'csspec_py2.log'
        spec['chatter']   = 3

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
        self.test_value(spectrum.nrows(), bins, 'Check for %d rows in spectrum' % bins)
        self.test_value(spectrum.ncols(), 8, 'Check for 8 columns in spectrum')

        # Return
        return
