#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csobsinfo script.
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
# Test class for csobsinfo script #
# =============================== #
class Test(test):
    """
    Test class for csobsinfo script

    This test class makes unit tests for the csobsinfo script by using it
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
        self._inobs = self._datadir + '/obs_unbinned.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csobsinfo')

        # Append tests
        self.append(self._test_cmd, 'Test csobsinfo on command line')
        self.append(self._test_python, 'Test csobsinfo from Python')

        # Return
        return

    # Test csobsinfo on command line
    def _test_cmd(self):
        """
        Test csobsinfo on the command line
        """
        # Set script name
        csobsinfo = self._script('csobsinfo')

        # Setup csobsinfo command
        cmd = csobsinfo+' inobs="'+self._inobs+'"'+ \
                        ' outds9file="csobsinfo_cmd1.reg"'+ \
                        ' offset=yes ra=83.63 dec=22.01'+ \
                        ' logfile="csobsinfo_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check DS9 file
        self._check_ds9_file('csobsinfo_cmd1.reg')

        # Setup csobsinfo command
        cmd = csobsinfo+' inobs="obs_definition_that_does_not_exist.xml"'+ \
                        ' outds9file="csobsinfo_cmd2.reg"'+ \
                        ' offset=yes ra=83.63 dec=22.01'+ \
                        ' logfile="csobsinfo_cmd2.log"'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csobsinfo from Python
    def _test_python(self):
        """
        Test csobsinfo from Python
        """
        # Set-up csobsinfo
        obsinfo = cscripts.csobsinfo()
        obsinfo['inobs']      = self._inobs
        obsinfo['outds9file'] = 'csobsinfo_py1.reg'
        obsinfo['offset']     = True
        obsinfo['ra']         = 83.63
        obsinfo['dec']        = 22.01
        obsinfo['logfile']    = 'csobsinfo_py1.log'
        obsinfo['chatter']    = 2

        # Run csobsinfo script and save XML file
        obsinfo.logFileOpen()   # Make sure we get a log file
        obsinfo.run()
        obsinfo.save()

        # Check DS9 file
        self._check_ds9_file('csobsinfo_py1.reg')

        # Set-up csobsinfo
        obsinfo = cscripts.csobsinfo()
        obsinfo['inobs']      = self._inobs
        obsinfo['outds9file'] = 'csobsinfo_py2.reg'
        obsinfo['offset']     = True
        obsinfo['ra']         = 83.63
        obsinfo['dec']        = 22.01
        obsinfo['logfile']    = 'csobsinfo_py2.log'
        obsinfo['chatter']    = 3

        # Execute csobsinfo script
        obsinfo.execute()

        # Check DS9 file
        self._check_ds9_file('csobsinfo_py2.reg')

        # Return
        return

    # Check DS9 file
    def _check_ds9_file(self, filename):
        """
        Check DS9 file
        """
        # Open file   
        f = open(filename,'r')

        # Expect "fk5" in first line
        line = f.readline().strip('\n')
        self.test_assert(line == 'fk5',
                         'Check for "fk5" in first line')

        # Expect "fk5" in first line
        line = f.readline().strip('\n')
        self.test_assert(line[0:6] == 'point(',
                         'Check for "point(" in second line')

        # Close file
        f.close()

        # Return
        return
