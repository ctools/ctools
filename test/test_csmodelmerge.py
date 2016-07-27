#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csmodelmerge script.
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


# ================================== #
# Test class for csmodelmerge script #
# ================================== #
class Test(test):
    """
    Test class for csmodelmerge script

    This test class makes unit tests for the csmodelmerge script by using it
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
        self._model1  = self._datadir + '/crab.xml'
        self._model2  = self._datadir + '/model_cube_background.xml'
        self._model3  = self._datadir + '/model_cube_background*.xml'
        self._model4  = '@' + self._datadir + '/models.txt'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('csmodelmerge')

        # Append tests
        self.append(self._test_cmd, 'Test csmodelmerge on command line')
        self.append(self._test_python, 'Test csmodelmerge from Python')

        # Return
        return

    # Test csmodelmerge on command line
    def _test_cmd(self):
        """
        Test csmodelmerge on the command line
        """
        # Set script name
        csmodelmerge = self._script('csmodelmerge')

        # Setup csmodelmerge command
        cmd = csmodelmerge+' inmodels="'+self._model1+' '+self._model2+'"'+ \
                           ' outmodel="csmodelmerge_cmd1.xml"'+ \
                           ' logfile="csmodelmerge_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check model file
        self._check_model_file('csmodelmerge_cmd1.xml', 3)

        # Setup csmodelmerge command
        cmd = csmodelmerge+' inmodels="model_that_does_not_exist.xml"'+ \
                           ' outmodel="csmodelmerge_cmd2.xml"'+ \
                           ' logfile="csmodelmerge_cmd2.log"'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csmodelmerge from Python
    def _test_python(self):
        """
        Test csmodelmerge from Python
        """
        # Set-up csmodelmerge for space-separated input
        modelmerge = cscripts.csmodelmerge()
        modelmerge['inmodels'] = self._model1+' '+self._model2
        modelmerge['outmodel'] = 'csmodelmerge_py1.xml'
        modelmerge['logfile']  = 'csmodelmerge_py1.log'
        modelmerge['chatter']  = 2

        # Run csmodelmerge script and save models
        modelmerge.logFileOpen()   # Make sure we get a log file
        modelmerge.run()
        modelmerge.save()

        # Check model file
        self._check_model_file('csmodelmerge_py1.xml', 3)

        # Set-up csmodelmerge for semi-colon separated input
        modelmerge = cscripts.csmodelmerge()
        modelmerge['inmodels'] = self._model1+';'+self._model2
        modelmerge['outmodel'] = 'csmodelmerge_py2.xml'
        modelmerge['logfile']  = 'csmodelmerge_py2.log'
        modelmerge['chatter']  = 3

        # Execute csmodelmerge script
        modelmerge.execute()

        # Check model file
        self._check_model_file('csmodelmerge_py2.xml', 3)

        # Set-up csmodelmerge for wildcard input
        modelmerge = cscripts.csmodelmerge()
        modelmerge['inmodels'] = self._model3
        modelmerge['outmodel'] = 'csmodelmerge_py3.xml'
        modelmerge['logfile']  = 'csmodelmerge_py3.log'
        modelmerge['chatter']  = 4

        # Execute csmodelmerge script
        modelmerge.execute()

        # Check model file
        self._check_model_file('csmodelmerge_py3.xml', 2)

        # Set-up csmodelmerge for ASCII file
        modelmerge = cscripts.csmodelmerge()
        modelmerge['inmodels'] = self._model4
        modelmerge['outmodel'] = 'csmodelmerge_py4.xml'
        modelmerge['logfile']  = 'csmodelmerge_py4.log'
        modelmerge['chatter']  = 4

        # Execute csmodelmerge script
        modelmerge.execute()

        # Check model file
        self._check_model_file('csmodelmerge_py4.xml', 2)

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
        self.test_value(models.size(), number,
                        'Check for '+str(number)+' models in XML file')
        
        # Return
        return
