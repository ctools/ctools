#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cssrcdetect script
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
import os
import gammalib
import cscripts
from testing import test


# ================================= #
# Test class for cssrcdetect script #
# ================================= #
class Test(test):
    """
    Test class for cssrcdetect script

    This test class makes unit tests for the cssrcdetect script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        test.__init__(self)

        # Set members
        self._map = self._datadir + '/skymap.fits'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('cssrcdetect')

        # Append tests
        self.append(self._test_cmd, 'Test cssrcdetect on command line')
        self.append(self._test_python, 'Test cssrcdetect from Python')

        # Return
        return

    # Test cssrcdetect on command line
    def _test_cmd(self):
        """
        Test cssrcdetect on the command line
        """
        # Set script name
        cssrcdetect = self._script('cssrcdetect')

        # Setup cssrcdetect command
        cmd = cssrcdetect+' inmap="'+self._map+'"'+ \
                          ' outmodel="cssrcdetect_cmd1.xml"'+ \
                          ' outds9file="cssrcdetect_cmd1.reg"'+ \
                          ' srcmodel="POINT" bkgmodel="NONE" threshold=10.0'+ \
                          ' logfile="cssrcdetect_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check model file
        self._check_model_file('cssrcdetect_cmd1.xml', 2)

        # Setup cssrcdetect command
        cmd = cssrcdetect+' inmap="skymap_that_does_not_exist"'+ \
                          ' outmodel="cssrcdetect_cmd1.xml"'+ \
                          ' outds9file="cssrcdetect_cmd1.reg"'+ \
                          ' srcmodel="POINT" bkgmodel="NONE" threshold=10.0'+ \
                          ' logfile="cssrcdetect_cmd2.log" chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test cssrcdetect from Python
    def _test_python(self):
        """
        Test cssrcdetect from Python
        """
        # Allocate empty cssrcdetect script
        srcdetect = cscripts.cssrcdetect()

        # Check that saving does not nothing
        srcdetect['outmodel']   = 'cssrcdetect_py0.xml'
        srcdetect['outds9file'] = 'cssrcdetect_py0.reg'
        srcdetect['logfile']    = 'cssrcdetect_py0.log'
        srcdetect.logFileOpen()
        srcdetect.save()

        # Check than an empty model definition file has been created
        self._check_model_file('cssrcdetect_py0.xml', 0)

        # Check that clearing does not lead to an exception or segfault
        srcdetect.clear()

        # Now set cssrcdetect parameters
        srcdetect['inmap']      = self._map
        srcdetect['outmodel']   = 'cssrcdetect_py1.xml'
        srcdetect['outds9file'] = 'cssrcdetect_py1.reg'
        srcdetect['srcmodel']   = 'POINT'
        srcdetect['bkgmodel']   = 'NONE'
        srcdetect['threshold']  = 10.0
        srcdetect['logfile']    = 'cssrcdetect_py1.log'
        srcdetect['chatter']    = 2

        # Run cssrcdetect script and save models
        srcdetect.logFileOpen()   # Make sure we get a log file
        srcdetect.run()
        srcdetect.save()

        # Check model file
        self._check_model_file('cssrcdetect_py1.xml', 2)

        # Execute cssrcdetect script again, now with a higher chatter level
        # than before
        srcdetect['outmodel']   = 'cssrcdetect_py2.xml'
        srcdetect['outds9file'] = 'cssrcdetect_py2.reg'
        srcdetect['bkgmodel']   = 'IRF'
        srcdetect['logfile']    = 'cssrcdetect_py2.log'
        srcdetect['chatter']    = 3
        srcdetect.logFileOpen()  # Needed to get a new log file
        srcdetect.execute()

        # Check model file
        self._check_model_file('cssrcdetect_py2.xml', 3)

        # Return
        return

    # Check model file
    def _check_model_file(self, filename, number):
        """
        Check model file
        """
        # Open models
        models = gammalib.GModels(filename)

        # Check number of models
        self.test_value(models.size(), number, 'Check for '+str(number)+
                        ' models in XML file')
        
        # Return
        return
