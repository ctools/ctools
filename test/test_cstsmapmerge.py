#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cstsmapmerge script.
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


# ================================== #
# Test class for cstsmapmerge script #
# ================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for cstsmapmerge script.

    This test class makes unit tests for the cstsmapmerge script by using it
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
        self.name('cstsmapmerge')

        # Append tests
        self.append(self._test_cmd, 'Test cstsmapmerge on command line')
        self.append(self._test_python, 'Test cstsmapmerge from Python')

        # Return
        return

    # Test cstsmapmerge on command line
    def _test_cmd(self):
        """
        Test cstsmapmerge on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile('README'):
            cstsmapmerge = '../cscripts/cstsmapmerge.py'
        else:
            cstsmapmerge = 'cstsmapmerge'

        # Setup cstsmapmerge command
        cmd = cstsmapmerge+' inmaps="data/tsmap_0.fits;data/tsmap_1.fits"'+ \
                           ' outmap="cstsmapmerge_cmd1.fits"'+ \
                           ' logfile="cstsmapmerge_cmd1.log" chatter=1'

        # Execute cstsmapmerge, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         'Successful cstsmapmerge execution on command line')

        # Check result file
        self._check_result_file('cstsmapmerge_cmd1.fits')

        # Setup cstsmapmerge command
        cmd = cstsmapmerge+' inmaps="map_that_does_not_exist.fits"'+ \
                           ' outmap="cstsmapmerge_cmd2.fits"'+ \
                           ' logfile="cstsmapmerge_cmd2.log" chatter=2'

        # Execute cstsmapmerge, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         'Failure of cstsmapmerge execution on command line')

        # Return
        return

    # Test cstsmapmerge from Python
    def _test_python(self):
        """
        Test cstsmapmerge from Python.
        """
        # Set-up semi-colon separated cstsmapmerge
        merge = cscripts.cstsmapmerge()
        merge['inmaps']  = 'data/tsmap_0.fits;data/tsmap_1.fits'
        merge['outmap']  = 'cstsmapmerge_py1.fits'
        merge['logfile'] = 'cstsmapmerge_py1.log'
        merge['chatter'] = 2

        # Run cstsmapmerge script
        merge.logFileOpen()   # Make sure we get a log file
        merge.run()
        merge.save()

        # Check pull distribution file
        self._check_result_file('cstsmapmerge_py1.fits')

        # Set-up space-separated cstsmapmerge
        merge = cscripts.cstsmapmerge()
        merge['inmaps']  = 'data/tsmap_0.fits data/tsmap_1.fits'
        merge['outmap']  = 'cstsmapmerge_py2.fits'
        merge['logfile'] = 'cstsmapmerge_py2.log'
        merge['chatter'] = 3

        # Run cstsmapmerge script
        merge.execute()

        # Check pull distribution file
        self._check_result_file('cstsmapmerge_py2.fits')

        # Set-up wildcard string cstsmapmerge
        merge = cscripts.cstsmapmerge()
        merge['inmaps']  = 'data/tsmap_*.fits'
        merge['outmap']  = 'cstsmapmerge_py3.fits'
        merge['logfile'] = 'cstsmapmerge_py3.log'
        merge['chatter'] = 4

        # Run cstsmapmerge script
        merge.execute()

        # Check pull distribution file
        self._check_result_file('cstsmapmerge_py3.fits')

        # Set-up ASCII file cstsmapmerge
        merge = cscripts.cstsmapmerge()
        merge['inmaps']  = '@data/tsmaps.txt'
        merge['outmap']  = 'cstsmapmerge_py4.fits'
        merge['logfile'] = 'cstsmapmerge_py4.log'
        merge['chatter'] = 4

        # Run cstsmapmerge script
        merge.execute()

        # Check pull distribution file
        self._check_result_file('cstsmapmerge_py4.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file.
        """
        # Open result file
        fits = gammalib.GFits(filename)

        # Get HDUs
        ts        = fits['Primary']
        prefactor = fits['Prefactor']
        index     = fits['Index']

        # Check dimensions
        self.test_value(ts.naxis(), 2, 'Check for 2 dimensions')
        self.test_value(ts.naxes(0), 2, 'Check for 2 pixels in X')
        self.test_value(ts.naxes(1), 1, 'Check for 1 pixel in Y')
        self.test_value(prefactor.naxis(), 2, 'Check for 2 dimensions')
        self.test_value(prefactor.naxes(0), 2, 'Check for 2 pixels in X')
        self.test_value(prefactor.naxes(1), 1, 'Check for 1 pixel in Y')
        self.test_value(index.naxis(), 2, 'Check for 2 dimensions')
        self.test_value(index.naxes(0), 2, 'Check for 2 pixels in X')
        self.test_value(index.naxes(1), 1, 'Check for 1 pixel in Y')

        # Return
        return
