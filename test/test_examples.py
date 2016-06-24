#! /usr/bin/env python
# ==========================================================================
# Performs unit tests for the ctools example scripts
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
import sys
import gammalib
import ctools
import cscripts


# =============================================== #
# Test class for example executables unit testing #
# =============================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib example executables
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

    # Execute Python script
    def _execute_python(self, name):
        """
        Execute Python script
        
        Parameters
        ----------
        name : str
            Python script name with .py extension.
        """
        # Setup command
        cmd = '../examples/' + name + '.py'

        # Execute Python script, make sure we catch any exception
        try:
            rc = os.system(cmd+' > example_'+name+'.log 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0, 'Check "'+name+'" script from command line')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('examples')

        # Append tests
        self.append(self.test_generate_prod3_irfs, 'Test generate_prod3_irfs')
        self.append(self.test_pipeline_binned_disk, 'Test pipeline_binned_disk')
        self.append(self.test_pipeline_binned_mem, 'Test pipeline_binned_mem')
        self.append(self.test_pipeline_stacked_disk, 'Test pipeline_stacked_disk')
        self.append(self.test_pipeline_stacked_mem, 'Test pipeline_stacked_mem')
        self.append(self.test_pipeline_unbinned_disk, 'Test pipeline_unbinned_disk')
        self.append(self.test_pipeline_unbinned_mem, 'Test pipeline_unbinned_mem')

        # Return
        return

    # Test generate_prod3_irfs
    def test_generate_prod3_irfs(self):
        """
        Test generate_prod3_irfs
        """
        # Execute script
        self._execute_python('generate_prod3_irfs')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_binned_disk
    def test_pipeline_binned_disk(self):
        """
        Test pipeline_binned_disk
        """
        # Execute script
        self._execute_python('pipeline_binned_disk')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_binned_mem
    def test_pipeline_binned_mem(self):
        """
        Test pipeline_binned_mem
        """
        # Execute script
        self._execute_python('pipeline_binned_mem')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_stacked_disk
    def test_pipeline_stacked_disk(self):
        """
        Test pipeline_stacked_disk
        """
        # Execute script
        self._execute_python('pipeline_stacked_disk')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_stacked_mem
    def test_pipeline_stacked_mem(self):
        """
        Test pipeline_stacked_mem
        """
        # Execute script
        self._execute_python('pipeline_stacked_mem')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_unbinned_disk
    def test_pipeline_unbinned_disk(self):
        """
        Test pipeline_unbinned_disk
        """
        # Execute script
        self._execute_python('pipeline_unbinned_disk')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_unbinned_mem
    def test_pipeline_unbinned_mem(self):
        """
        Test pipeline_unbinned_mem
        """
        # Execute script
        self._execute_python('pipeline_unbinned_mem')

        #TODO: Do any testing

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Set calibration database
    os.environ['CALDB'] = '../caldb'

    # Set PFILES environment variable
    try:
        os.mkdir('pfiles')
    except:
        pass
    os.environ['PFILES'] = 'pfiles'

    # Copy over pfiles
    os.system('cp -r ../src/*/*.par pfiles/')
    os.system('cp -r ../cscripts/*.par pfiles/')

    # Allocate test suites
    suites = gammalib.GTestSuites('Examples testing')

    # Allocate test suite, setup tests and append them to the container
    suite = Test()
    suite.set()
    suites.append(suite)

    # Run test suite
    success = suites.run()

    # Save test results
    suites.save('reports/examples.xml')

    # Set return code
    if success:
        rc = 0
    else:
        rc = 1

    # Exit with return code
    sys.exit(rc)
