#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csmodelselect script
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


# =================================== #
# Test class for csmodelselect script #
# =================================== #
class Test(test):
    """
    Test class for csmodelselect script

    This test class makes unit tests for the csmodelselect script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor.
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
        self.name('csmodelselect')

        # Append tests
        self.append(self._test_cmd, 'Test csmodelselect on command line')
        self.append(self._test_python, 'Test csmodelselect from Python')

        # Return
        return

    # Test csmodelselect on command line
    def _test_cmd(self):
        """
        Test csmodelselect on the command line
        """
        # Set script name
        csmodelselect = self._script('csmodelselect')

        # Setup csmodelselect command
        cmd = csmodelselect+' inobs="'+self._events+'"'+ \
                            ' inmodel="'+self._model+'"'+ \
                            ' outmodel="csmodelselect_cmd1.xml"'+ \
                            ' logfile="csmodelselect_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check model file
        self._check_model_file('csmodelselect_cmd1.xml', 2)

        # Setup csmodelselect command
        cmd = csmodelselect+' inobs="'+self._events+'"'+ \
                            ' inmodel="model_that_does_not_exist.xml"'+ \
                            ' outmodel="csmodelselect_cmd2.xml"'+ \
                            ' logfile="csmodelselect_cmd2.log"'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csmodelselect from Python
    def _test_python(self):
        """
        Test csmodelselect from Python
        """
        # Allocate empty csmodelselect script
        modelselect = cscripts.csmodelselect()

        # Check that saving does not nothing
        modelselect['outmodel'] = 'csmodelselect_py0.xml'
        modelselect['logfile']  = 'csmodelselect_py0.log'
        modelselect.logFileOpen()
        modelselect.save()

        # Check than an empty model definition file has been created
        self._check_model_file('csmodelselect_py0.xml', 0)

        # Check that clearing does not lead to an exception or segfault
        modelselect.clear()

        # Now set csmodelselect parameters
        modelselect['inobs']    = self._events
        modelselect['inmodel']  = self._model
        modelselect['outmodel'] = 'csmodelselect_py1.xml'
        modelselect['logfile']  = 'csmodelselect_py1.log'
        modelselect['chatter']  = 2

        # Run csmodelselect script and save models
        modelselect.logFileOpen()   # Make sure we get a log file
        modelselect.run()
        modelselect.save()

        # Check model file
        self._check_model_file('csmodelselect_py1.xml', 2)

        # Execute csmodelselect script again, now with a higher chatter level
        # than before
        modelselect['outmodel'] = 'csmodelselect_py2.xml'
        modelselect['logfile']  = 'csmodelselect_py2.log'
        modelselect['chatter']  = 3
        modelselect.logFileOpen()  # Needed to get a new log file
        modelselect.execute()

        # Check model file
        self._check_model_file('csmodelselect_py2.xml', 2)

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
