#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csworkflow script.
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


# ================================ #
# Test class for csworkflow script #
# ================================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for csworkflow script.

    This test class makes unit tests for the cslightcrv script by using it
    from the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
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

        # Kluge to set the command (installed version has no README file)
        if os.path.isfile('README'):
            csworkflow = '../cscripts/csworkflow.py'
        else:
            csworkflow = 'csworkflow'

        # Remove result file
        gammalib.GFilename('wf_crab_results.xml').remove()

        # Setup csworkflow command
        cmd = csworkflow+' inflow="data/workflow.xml"'+ \
                         ' logfile="csworkflow_cmd1.log" chatter=1'

        # Execute csworkflow, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         'Successful csworkflow execution on command line')

        # Check fit result
        self._check_fit_result('wf_crab_results.xml')

        # Return
        return

    # Test csworkflow from Python
    def _test_python(self):

        # Remove result file
        gammalib.GFilename('wf_crab_results.xml').remove()

        # Set-up csworkflow
        workflow = cscripts.csworkflow()
        workflow['inflow']  = 'data/workflow.xml'
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
        workflow['inflow']  = 'data/workflow.xml'
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

        # Load fit results
        models = gammalib.GModels(filename)

        # Check fit result values
        self.test_value(models['Crab'][2].value(),
                        5.73787e-16, 1.0e-5, 'Crab prefactor')
        self.test_value(models['Crab'][3].value(),
                        -2.47386, 1.0e-5, 'Crab index')
        self.test_value(models['CTABackgroundModel'][0].value(),
                        1.01414, 1.0e-5, 'Background model prefactor')
        self.test_value(models['CTABackgroundModel'][1].value(),
                        0.00242827, 1.0e-5, 'Background model index')

        # Return
        return
