#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cstsdist script.
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


# ============================== #
# Test class for cstsdist script #
# ============================== #
class Test(test):
    """
    Test class for cstsdist script

    This test class makes unit tests for the cstsdist script by using it
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
        self.name('cstsdist')

        # Append tests
        self.append(self._test_cmd, 'Test cstsdist on command line')
        self.append(self._test_python, 'Test cstsdist from Python')

        # Return
        return

    # Test cstsdist on command line
    def _test_cmd(self):
        """
        Test cstsdist on the command line
        """
        # Set script name
        cstsdist = self._script('cstsdist')

        # Setup cstsdist command
        cmd = cstsdist+' inmodel="'+self._model+'"'+ \
                       ' srcname="Crab"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                       ' ntrials=1 ra=83.63 dec=22.01 emin=0.1 emax=100.0'+ \
                       ' enumbins=0 tmax=1800.0 rad=5.0 npix=200 npix=200'+ \
                       ' binsz=0.05'+ \
                       ' outfile="cstsdist_cmd1.dat"'+ \
                       ' logfile="cstsdist_cmd1.log" chatter=2'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('cstsdist_cmd1.dat')

        # Setup cstsdist command
        cmd = cstsdist+' inmodel="model_that_does_not_exist.xml"'+ \
                       ' srcname="Crab"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                       ' ntrials=1 ra=83.63 dec=22.01 emin=0.1 emax=100.0'+ \
                       ' enumbins=0 tmax=1800.0 rad=5.0 npix=200 npix=200'+ \
                       ' binsz=0.05'+ \
                       ' outfile="cstsdist_cmd2.dat"'+ \
                       ' logfile="cstsdist_cmd2.log" chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test cstsdist from Python
    def _test_python(self):
        """
        Test cstsdist from Python
        """
        # Set-up cstsdist
        tsdist = cscripts.cstsdist()
        tsdist['inmodel']  = self._model
        tsdist['srcname']  = 'Crab'
        tsdist['caldb']    = self._caldb
        tsdist['irf']      = self._irf
        tsdist['ntrials']  = 1
        tsdist['ra']       = 83.63
        tsdist['dec']      = 22.01
        tsdist['emin']     = 0.1
        tsdist['emax']     = 100.0
        tsdist['enumbins'] = 0
        tsdist['tmax']     = 1800.0
        tsdist['rad']      = 5.0
        tsdist['outfile']  = 'cstsdist_py1.dat'
        tsdist['logfile']  = 'cstsdist_py1.log'
        tsdist['chatter']  = 3

        # Run cstsdist script
        tsdist.logFileOpen()   # Make sure we get a log file
        tsdist.run()
        tsdist.save()

        # Check pull distribution file
        self._check_result_file('cstsdist_py1.dat')

        # Set-up binned cstsdist in four observation wobble mode
        tsdist = cscripts.cstsdist()
        tsdist['inmodel']  = self._model
        tsdist['srcname']  = 'Crab'
        tsdist['caldb']    = self._caldb
        tsdist['irf']      = self._irf
        tsdist['ntrials']  = 1
        tsdist['ra']       = 83.63
        tsdist['dec']      = 22.01
        tsdist['coordsys'] = 'CEL'
        tsdist['proj']     = 'TAN'
        tsdist['emin']     = 0.1
        tsdist['emax']     = 100.0
        tsdist['enumbins'] = 10
        tsdist['tmax']     = 1800.0
        tsdist['rad']      = 5.0
        tsdist['npix']     = 50
        tsdist['binsz']    = 0.1
        tsdist['pattern']  = 'four'
        tsdist['outfile']  = 'cstsdist_py2.dat'
        tsdist['logfile']  = 'cstsdist_py2.log'
        tsdist['chatter']  = 4

        # Execute cstsdist script
        tsdist.execute()

        # Check pull distribution file
        self._check_result_file('cstsdist_py2.dat', ncols=9)

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, ncols=11):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Name of result file
        ncols : int, optional
            Required number of columns
        """
        # Open result file as CSV file
        results = gammalib.GCsv(filename, ',')

        # Check dimensions
        self.test_value(results.nrows(), 2,     'Check rows in TS file')
        self.test_value(results.ncols(), ncols, 'Check columns in TS file')

        # Return
        return
