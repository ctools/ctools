#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctulimit tool.
#
# Copyright (C) 2015-2021 Michael Mayer
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
import ctools
from testing import test


# ============================ #
# Test class for ctulimit tool #
# ============================ #
class Test(test):
    """
    Test class for ctulimit tool
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
        self.name('ctulimit')

        # Append tests
        self.append(self._test_cmd, 'Test ctulimit on command line')
        self.append(self._test_python, 'Test ctulimit from Python')

        # Return
        return

    # Test ctulimit on command line
    def _test_cmd(self):
        """
        Test ctulimit on the command lines
        """
        # Set tool name
        ctulimit = self._tool('ctulimit')

        # Setup ctulimit command
        cmd = ctulimit+' inobs="'+self._events+'"'+ \
                       ' inmodel="'+self._model+'" srcname="Crab"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                       ' tol=0.1 logfile="ctulimit_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Setup ctulimit command
        cmd = ctulimit+' inobs="event_file_that_does_not_exist.fits"'+ \
                       ' inmodel="'+self._model+'" srcname="Crab"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                       ' tol=0.1 logfile="ctulimit_cmd2.log" debug=yes chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctulimit command
        cmd = ctulimit+' inobs="'+self._events+'"'+ \
                       ' inmodel="'+self._model+'" srcname="Vela"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                       ' tol=0.1 logfile="ctulimit_cmd2.log" debug=yes chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid source name when executed from command line')

        # Check ctulimit --help
        self._check_help(ctulimit)

        # Return
        return

    # Test ctulimit from Python
    def _test_python(self):
        """
        Test ctulimit from Python
        """
        # Allocate ctulimit
        ulimit = ctools.ctulimit()

        # Check that empty ctulimit tool holds zero upper limits
        self.test_value(ulimit.diff_ulimit(), 0.0, 1.0e-21,
                        'Check differential upper limit')
        self.test_value(ulimit.flux_ulimit(), 0.0, 1.0e-16,
                        'Check upper limit on photon flux')
        self.test_value(ulimit.eflux_ulimit(), 0.0, 1.0e-16,
                        'Check upper limit on energy flux')

        # Check that saving does not nothing
        ulimit['logfile'] = 'ctulimit_py0.log'
        ulimit.logFileOpen()
        ulimit.save()

        # Check that clearing does not lead to an exception or segfault
        ulimit.clear()

        # Now set ctulimit parameters
        ulimit['inobs']   = self._events
        ulimit['inmodel'] = self._model
        ulimit['srcname'] = 'Crab'
        ulimit['caldb']   = self._caldb
        ulimit['irf']     = self._irf
        ulimit['tol']     = 0.1
        ulimit['logfile'] = 'ctulimit_py1.log'
        ulimit['chatter'] = 2

        # Run ctulimit tool
        ulimit.logFileOpen()   # Make sure we get a log file
        ulimit.run()
        ulimit.save()

        # Set reference value
        ref_diff  = 2.75359655874408e-17
        ref_flux  = 1.66942609234064e-11
        ref_eflux = 6.4588860064421e-11

        # Check results
        self.test_value(ulimit.diff_ulimit(), ref_diff, 1.0e-21,
                        'Check differential upper limit')
        self.test_value(ulimit.flux_ulimit(), ref_flux, 1.0e-16,
                        'Check upper limit on photon flux')
        self.test_value(ulimit.eflux_ulimit(), ref_eflux, 1.0e-16,
                        'Check upper limit on energy flux')

        # Check obs() method
        self.test_value(ulimit.obs().size(), 1,
                        'Check number of observations in container')

        # Check opt() method
        self.test_value(ulimit.opt().status(), 0, 'Check optimizer status')

        # Copy ctulimit tool
        cpy_ulimit = ulimit.copy()

        # Check results of copy
        self.test_value(cpy_ulimit.diff_ulimit(), ref_diff, 1.0e-21,
                        'Check differential upper limit')
        self.test_value(cpy_ulimit.flux_ulimit(), ref_flux, 1.0e-16,
                        'Check upper limit on photon flux')
        self.test_value(cpy_ulimit.eflux_ulimit(), ref_eflux, 1.0e-16,
                        'Check upper limit on energy flux')

        # Now clear copy of ctulimit tool
        cpy_ulimit.clear()

        # Check that empty ctulimit tool holds zero upper limits
        self.test_value(cpy_ulimit.diff_ulimit(), 0.0, 1.0e-21,
                        'Check differential upper limit')
        self.test_value(cpy_ulimit.flux_ulimit(), 0.0, 1.0e-16,
                        'Check upper limit on photon flux')
        self.test_value(cpy_ulimit.eflux_ulimit(), 0.0, 1.0e-16,
                        'Check upper limit on energy flux')

        # Run ctlike to get an initial log-likelihood solution
        like = ctools.ctlike()
        like['inobs']    = self._events
        like['inmodel']  = self._model
        like['caldb']    = self._caldb
        like['irf']      = self._irf
        like.run()

        # Now set ctulimit tool using the observation container from the
        # previous run. This should avoid the necessity to recompute the
        # maximum likelihood
        ulimit = ctools.ctulimit(like.obs())
        ulimit['srcname'] = 'Crab'
        ulimit['tol']     = 0.1
        ulimit['logfile'] = 'ctulimit_py2.log'
        ulimit['chatter'] = 3

        # Execute ctulimit tool
        ulimit.logFileOpen()  # Needed to get a new log file
        ulimit.execute()

        # Check results
        self.test_value(ulimit.diff_ulimit(), ref_diff, 1.0e-21,
                        'Check differential upper limit')
        self.test_value(ulimit.flux_ulimit(), ref_flux, 1.0e-16,
                        'Check upper limit on photon flux')
        self.test_value(ulimit.eflux_ulimit(), ref_eflux, 1.0e-16,
                        'Check upper limit on energy flux')

        # Test invalid model name
        ulimit['srcname'] = 'Weihnachtsstern'
        ulimit['logfile'] = 'ctulimit_py3.log'
        ulimit.logFileOpen()
        self.test_try('Test invalid model name')
        try:
            ulimit.execute()
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test specification of background model
        ulimit['srcname'] = 'Background'
        ulimit['logfile'] = 'ctulimit_py4.log'
        ulimit.logFileOpen()
        self.test_try('Test invalid model name')
        try:
            ulimit.execute()
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Test run with too few iterations
        ulimit['srcname']  = 'Crab'
        ulimit['max_iter'] = 1
        ulimit['logfile']  = 'ctulimit_py5.log'
        ulimit.logFileOpen()
        self.test_try('Test ctulimit with too few iterations')
        try:
            ulimit.execute()
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Return
        return
