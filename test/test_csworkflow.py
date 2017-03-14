#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csworkflow script.
#
# Copyright (C) 2016-2017 Juergen Knoedlseder
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


# ================================ #
# Test class for csworkflow script #
# ================================ #
class Test(test):
    """
    Test class for csworkflow script

    This test class makes unit tests for the cslightcrv script by using it
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
        self._workflow = self._datadir + '/workflow.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csworkflow')

        # Append tests
        self.append(self._test_cmd, 'Test csworkflow on command line')
        self.append(self._test_python, 'Test csworkflow from Python')

        # Return
        return

    # Test csworkflow on command line
    def _test_cmd(self):
        """
        Test csworkflow on the command line
        """
        # Set script name
        csworkflow = self._script('csworkflow')

        # Remove result file
        gammalib.GFilename('wf_crab_results.xml').remove()

        # Setup csworkflow command
        cmd = csworkflow+' inflow="'+self._workflow+'"'+ \
                         ' logfile="csworkflow_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check fit result
        self._check_fit_result('wf_crab_results.xml')

        # Return
        return

    # Test csworkflow from Python
    def _test_python(self):
        """
        Test csworkflow from Python
        """
        # Remove result file
        gammalib.GFilename('wf_crab_results.xml').remove()

        # Set-up csworkflow
        workflow = cscripts.csworkflow()
        workflow['inflow']  = self._workflow
        workflow['logfile'] = 'csworkflow_py1.log'
        workflow['chatter'] = 2

        # Run script
        workflow.logFileOpen()   # Make sure we get a log file
        workflow.run()

        # Check fit result
        self._check_fit_result('wf_crab_results.xml')

        # Remove result file
        gammalib.GFilename('wf_crab_results.xml').remove()

        # Set-up csworkflow
        workflow = cscripts.csworkflow()
        workflow['inflow']  = self._workflow
        workflow['logfile'] = 'csworkflow_py2.log'
        workflow['chatter'] = 3

        # Execute script
        workflow.execute()

        # Check fit result
        self._check_fit_result('wf_crab_results.xml')

        # Return
        return

    # Check fit result XML values
    def _check_fit_result(self, filename):
        """
        Check result file
        """
        # Load fit results
        models = gammalib.GModels(filename)

        # Check fit result values
        self.test_value(models['Crab'][2].value(),
                        5.71055e-16, 1.0e-5, 'Crab prefactor')
        self.test_value(models['Crab'][3].value(),
                        -2.463392, 1.0e-5, 'Crab index')
        self.test_value(models['Background'][0].value(),
                        3.008828, 1.0e-5, 'Background model sigma')
        self.test_value(models['Background'][1].value(),
                        6.27131e-05, 1.0e-5, 'Background model prefactor')
        self.test_value(models['Background'][2].value(),
                        -1.833184, 1.0e-5, 'Background model index')

        # Return
        return
