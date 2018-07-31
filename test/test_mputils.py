#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the mputils module
#
# Copyright (C) 2018 Juergen Knoedlseder
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
import ctools
from cscripts import mputils
from testing import test


# ================== #
# mputils test class #
# ================== #
class mputils_test(ctools.cscript):
    """
    Test class for the mputils module
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the base class constructor
        self._init_cscript(self.__class__.__name__, ctools.__version__, argv)

        # Return
        return

    # Bin computation
    def _compute_bin(self, i):
        """
        Bin computation
        """
        # Return bin index
        return i


# ============================= #
# Test class for mputils module #
# ============================= #
class Test(test):
    """
    Test class for mputils module
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Create parameter file for mputils_test application
        pars = gammalib.GApplicationPars()
        pars.append_standard()
        pars.save('mputils_test.par')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('mputils')

        # Append tests
        self.append(self._test_nthreads, 'Test mputils.nthreads() function')
        self.append(self._test_process, 'Test mputils.process() function')
        self.append(self._test_mpfunc, 'Test mputils.mpfunc() function')

        # Return
        return

    # Setup mputils test class
    def _setup_mputils_test(self):
        """
        Setup mputils test class
        """
        # Get test class
        cls = mputils_test()

        # Recover User parameters
        pars = cls.pars()

        # Create a copy of the User parameters
        cpy_pars = pars.copy()

        # Return
        return (cls, pars, cpy_pars)

    # Test nthreads() function
    def _test_nthreads(self):
        """
        Test nthreads() function
        """
        # Set reference number of CPU's
        ncpus = 1
        try:
            from multiprocessing import cpu_count
            ncpus = cpu_count()
        except:
            pass

        # Check nthreads=0. In this case the nthreads() function should return
        # the number of CPU's that are available, provided that the
        # multiprocessing module is available. If the module is not available
        # the function returns 1.
        cls, pars, cpy_pars = self._setup_mputils_test()
        pars.append(gammalib.GApplicationPar('nthreads','i','h','0','','',''))
        nthreads = mputils.nthreads(cls)
        self.test_value(nthreads, ncpus, 'Check nthreads=0')
        cls.pars(cpy_pars)

        # Check nthreads=739. In this case the nthreads() function should
        # return 739, provided that the multiprocessing module is available.
        # If the module is not available the function returns 1.
        cls, pars, cpy_pars = self._setup_mputils_test()
        pars.append(gammalib.GApplicationPar('nthreads','i','h','739','','',''))
        nthreads = mputils.nthreads(cls)
        self.test_value(nthreads, 739, 'Check nthreads=739')
        cls.pars(cpy_pars)

        # Return
        return

    # Test process() function
    def _test_process(self):
        """
        Test process() function
        """
        # Test only if Pool module is available
        try:
            # Import Pool module
            from multiprocessing import Pool

            # Get test class
            cls = mputils_test()

            # Setup arguments for mpfunc() computation
            args = [(cls, '_compute_bin', i) for i in range(10)]

            # Multiprocessing
            poolresults = mputils.process(2, mputils.mpfunc, args)

            # Loop over arguments
            for i in range(10):
                result = poolresults[i][0]
                self.test_value(result, i, 'Check process for argument %d' % i)

        # ... otherwise if Pool is not available then skip this test
        except:
            pass

        # Return
        return

    # Test mpfunc() function
    def _test_mpfunc(self):
        """
        Test mpfunc() function
        """
        # Get test class
        cls = mputils_test()

        # Setup arguments for mpfunc() computation
        args = [(cls, '_compute_bin', i) for i in range(10)]

        # Loop over arguments
        for arg in args:
            result, _ = mputils.mpfunc(arg)
            self.test_value(result, arg[2], 'Check mpfunc for argument %d' % arg[2])

        # Return
        return
