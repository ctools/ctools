#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the comobsselect script.
#
# Copyright (C) 2021 Juergen Knoedlseder
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
import comscripts
from testing import test


# ================================== #
# Test class for comobsselect script #
# ================================== #
class Test(test):
    """
    Test class for comobsselect script

    This test class makes unit tests for the comobsselect script by using it
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
        Set all test functions.
        """
        # Set test name
        self.name('comobsselect')

        # Append tests
        self.append(self._test_cmd, 'Test comobsselect on command line')
        self.append(self._test_python, 'Test comobsselect from Python')

        # Return
        return

    # Test command line use
    def _test_cmd(self):
        """
        Test command line use.
        """
        # Set script name
        comobsselect = self._comscript('comobsselect')

        # Setup command
        cmd = comobsselect+' pntselect="CIRCLE" coordsys="GAL" '+ \
                           'glon="184.5597" glat="-5.7892" rad="30.0" '+ \
                           'tmin="1991-05-01T00:00:00" '+ \
                           'tmax="1991-06-01T00:00:00" '+ \
                           'outobs="comobsselect_cmd1.xml" '+ \
                           'logfile="comobsselect_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result
        self._check_result('comobsselect_cmd1.xml')

        # Check --help option
        self._check_help(comobsselect)

        # Return
        return

    # Test from Python
    def _test_python(self):
        """
        Test from Python
        """
        # Set-up script
        select = comscripts.comobsselect()
        select['pntselect'] = 'CIRCLE'
        select['coordsys']  = 'GAL'
        select['glon']      = 184.5597
        select['glat']      = -5.7892
        select['rad']       = 30.0
        select['tmin']      = '1991-05-01T00:00:00'
        select['tmax']      = '1991-06-01T00:00:00'
        select['outobs']    = 'comobsselect_py1.xml'
        select['logfile']   = 'comobsselect_py1.log'
        select['chatter']   = 2

        # Run script and save result
        select.logFileOpen()   # Make sure we get a log file
        select.run()
        select.save()

        # Check result
        self._check_result('comobsselect_py1.xml')

        # Return
        return

    # Check result file
    def _check_result(self, filename, nobs=1):
        """
        Check result file
        """
        # Load observation definition file
        obs = gammalib.GObservations(filename)

        # Check that there is one observation
        self.test_value(obs.size(), nobs, 'Check for %d observations' % nobs)

        # Return
        return
