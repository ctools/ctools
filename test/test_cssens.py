#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cssens script
#
# Copyright (C) 2016-2021 Juergen Knoedlseder
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
# Test class for cssens script #
# ============================ #
class Test(test):
    """
    Test class for cssens script

    This test class makes unit tests for the cssens script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set test data
        self._model = self._datadir + '/model_crab_radialacceptance.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('cssens')

        # Append tests
        self.append(self._test_cmd, 'Test cssens on command line')
        self.append(self._test_python, 'Test cssens from Python')
        self.append(self._test_pickeling, 'Test cssens pickeling')

        # Return
        return

    # Test cssens on command line
    def _test_cmd(self):
        """
        Test cssens on the command line
        """
        # Set script name
        cssens = self._script('cssens')

        # Setup cssens command
        cmd = cssens+' inmodel="'+self._model+'"'+ \
                     ' srcname="Crab"'+ \
                     ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                     ' duration=180.0 rad=3.0 emin=1.0 emax=10.0 bins=1'+ \
                     ' outfile="cssens_cmd1.fits"'+ \
                     ' logfile="cssens_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('cssens_cmd1.fits', 1)

        # Setup cssens command
        cmd = cssens+' inmodel="model_that_does_not_exist.xml"'+ \
                     ' srcname="Crab"'+ \
                     ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                     ' duration=180.0 rad=3.0 emin=1.0 emax=10.0 bins=1'+ \
                     ' outfile="cssens_cmd2.fits"'+ \
                     ' logfile="cssens_cmd2.log" debug=yes chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check cssens --help
        self._check_help(cssens)

        # Return
        return

    # Test cssens from Python
    def _test_python(self):
        """
        Test cssens from Python
        """
        # Set-up cssens for differential sensitivity computation with
        # multiprocessing
        sens = cscripts.cssens()
        sens['inobs']    = 'NONE'
        sens['inmodel']  = self._model
        sens['srcname']  = 'Crab'
        sens['caldb']    = self._caldb
        sens['irf']      = self._irf
        sens['duration'] = 180.0
        sens['rad']      = 3.0
        sens['emin']     = 1.0
        sens['emax']     = 10.0
        sens['bins']     = 2
        sens['outfile']  = 'cssens_py1.fits'
        sens['logfile']  = 'cssens_py1.log'
        sens['chatter']  = 4

        # Run cssens script
        sens.logFileOpen()   # Make sure we get a log file
        sens.run()
        sens.save()

        # Check sensitivity file
        self._check_result_file('cssens_py1.fits', 2)

        # Set-up cssens for integral sensitivity computation without
        # multiprocessing
        sens = cscripts.cssens()
        sens['inobs']    = 'NONE'
        sens['inmodel']  = self._model
        sens['srcname']  = 'Crab'
        sens['caldb']    = self._caldb
        sens['irf']      = self._irf
        sens['duration'] = 180.0
        sens['rad']      = 3.0
        sens['emin']     = 1.0
        sens['emax']     = 10.0
        sens['bins']     = 1
        sens['type']     = 'Integral'
        sens['outfile']  = 'cssens_py2.fits'
        sens['logfile']  = 'cssens_py2.log'
        sens['chatter']  = 4
        sens['nthreads'] = 1

        # Execute cssens script
        sens.execute()

        # Check sensitivity file
        self._check_result_file('cssens_py2.fits', 1)

        # Return
        return

    # Test cssens pickeling
    def _test_pickeling(self):
        """
        Test cssens pickeling
        """
        # Perform pickeling tests of empty class
        self._pickeling(cscripts.cssens())

        # Set-up unbinned cssens
        sens = cscripts.cssens()
        sens['inobs']    = 'NONE'
        sens['inmodel']  = self._model
        sens['srcname']  = 'Crab'
        sens['caldb']    = self._caldb
        sens['irf']      = self._irf
        sens['duration'] = 180.0
        sens['rad']      = 3.0
        sens['emin']     = 1.0
        sens['emax']     = 10.0
        sens['bins']     = 2
        sens['outfile']  = 'cssens_py1_pickle.fits'
        sens['logfile']  = 'cssens_py1_pickle.log'
        sens['chatter']  = 4

        # Perform pickeling tests of filled class
        obj = self._pickeling(sens)

        # Run csspec script and save light curve
        obj.logFileOpen()   # Make sure we get a log file
        obj.run()
        obj.save()

        # Check result file
        self._check_result_file('cssens_py1_pickle.fits', 2)

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nrows):
        """
        Check result file
        """
        # Open result file as FITS file
        fits  = gammalib.GFits(filename)
        table = fits.table('SENSITIVITY')

        # Check dimensions
        self.test_value(table.nrows(), nrows,
             'Check number of rows in sensitivity file')
        self.test_value(table.ncols(), 10,
             'Check number of columns in sensitivity file')

        # Return
        return
