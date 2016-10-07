#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csviscube script
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


# =============================== #
# Test class for csviscube script #
# =============================== #
class Test(test):
    """
    Test class for csviscube script
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
        self.name('csviscube')

        # Append tests
        self.append(self._test_cmd, 'Test csviscube on command line')
        self.append(self._test_python, 'Test csviscube from Python')

        # Return
        return

    # Test csviscube on command line
    def _test_cmd(self):
        """
        Test csviscube on the command line
        """
        # Set script name
        csviscube = self._script('csviscube')

        # Setup csviscube command
        cmd = csviscube+' tmin=0.0 tmax=31557600.0 geolon=79.4041'+\
                        ' geolat=-24.6272 binsz=3.0 outfile=csviscube_cmd1.fits'+\
                        ' logfile="csviscube_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Return
        return

    # Test csviscube from Python
    def _test_python(self):
        """
        Test csviscube from Python
        """
        # Set-up csviscube
        viscube = cscripts.csviscube()
        viscube['tmin']    =  0.0
        viscube['tmax']    =  31557600.0
        viscube['geolon']  =  79.4041
        viscube['geolat']  = -24.6272
        viscube['binsz']   = 3.0
        viscube['outfile'] = 'csviscube_py1.fits'
        viscube['logfile'] = 'csviscube_py1.log'
        viscube['chatter'] = 2
        viscube['publish'] = True

        # Run script
        viscube.logFileOpen()   # Make sure we get a log file
        viscube.run()

        # Return
        return
