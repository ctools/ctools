#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctcubemask tool.
#
# Copyright (C) 2014-2016 Juergen Knoedlseder
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
import ctools
from testing import test


# ============================== #
# Test class for ctcubemask tool #
# ============================== #
class Test(test):
    """
    Test class for ctcubemask tool
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set members
        self._exclusion = self._datadir + '/exclusion.reg'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('ctcubemask')

        # Append tests
        self.append(self._test_cmd, 'Test ctcubemask on command line')
        self.append(self._test_python, 'Test ctcubemask from Python')

        # Return
        return

    # Test ctcubemask on command line
    def _test_cmd(self):
        """
        Test ctcubemask on the command line
        """
        # Set tool name
        ctcubemask = self._tool('ctcubemask')

        # Setup ctcubemask command
        cmd = ctcubemask+' inobs="'+self._cntcube+'"'+ \
                         ' regfile="'+self._exclusion+'"'+ \
                         ' outcube="ctcubemask_cmd1.fits"'+ \
                         ' ra=83.63 dec=22.01 rad=2.0'+ \
                         ' emin=0.1 emax=100.0'+ \
                         ' logfile="ctcubemask_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctcubemask_cmd1.fits')

        # Setup ctcubemask command
        cmd = ctcubemask+' inobs="counts_cube_that_does_not_exist.fits"'+ \
                         ' regfile="'+self._exclusion+'"'+ \
                         ' outcube="ctcubemask_cmd2.fits"'+ \
                         ' ra=83.63 dec=22.01 rad=2.0'+ \
                         ' emin=0.1 emax=100.0'+ \
                         ' logfile="ctcubemask_cmd2.log" chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctcubemask --help option
        cmd = ctcubemask+' --help'

        # Check if execution was successful in case that the CTOOLS
        # environment variable was set or failed otherwise
        if 'CTOOLS' in os.environ:
            self.test_value(self._execute(cmd), 0,
                 'Check successful execution with --help option')
        else:
            self.test_assert(self._execute(cmd) != 0,
                 'Check execution failure with --help option')

        # Return
        return

    # Test ctcubemask from Python
    def _test_python(self):
        """
        Test ctcubemask from Python
        """
        # Allocate ctcubemask
        mask = ctools.ctcubemask()

        # Check that empty ctcubemask tool holds an empty observation
        self._check_obs(mask.obs(), nobs=0)

        # Check that saving saves an empty counts cube
        mask['outcube'] = 'ctcubemask_py0.fits'
        mask['logfile'] = 'ctcubemask_py0.log'
        mask.logFileOpen()
        mask.save()
        self.test_assert(not os.path.isfile('ctcubemask_py0.fits'),
             'Check that no counts cube has been created')

        # Check that clearing does not lead to an exception or segfault
        mask.clear()

        # Now set ctcubemask parameters
        mask['inobs']   = self._cntcube
        mask['regfile'] = self._exclusion
        mask['ra']      = 83.63
        mask['dec']     = 22.01
        mask['rad']     = 2.0
        mask['emin']    = 0.1
        mask['emax']    = 100.0
        mask['outcube'] = 'ctcubemask_py1.fits'
        mask['logfile'] = 'ctcubemask_py1.log'
        mask['chatter'] = 2

        # Run ctcubemask tool
        mask.logFileOpen()   # Make sure we get a log file
        mask.run()
        mask.save()

        # Check result file
        self._check_result_file('ctcubemask_py1.fits')

        # Set-up ctcubemask without exclusion regions but tighter energy
        # selection
        mask = ctools.ctcubemask()
        mask['inobs']   = self._cntcube
        mask['regfile'] = 'NONE'
        mask['ra']      = 83.63
        mask['dec']     = 22.01
        mask['rad']     = 3.0
        mask['emin']    = 1.0
        mask['emax']    = 100.0
        mask['outcube'] = 'ctcubemask_py2.fits'
        mask['logfile'] = 'ctcubemask_py2.log'
        mask['chatter'] = 3
        mask['publish'] = True

        # Execute ctcubemask tool
        mask.logFileOpen()   # Make sure we get a log file
        mask.execute()

        # Check result file
        self._check_result_file('ctcubemask_py2.fits', events=573)

        # Copy ctcubemask tool
        cpy_mask = mask.copy()

        # Check that ctcubemask tool holds one observation
        self._check_obs(cpy_mask.obs())

        # Execute copy of ctcubemask tool again, now with a higher chatter
        # level than before
        cpy_mask['outcube'] = 'ctcubemask_py3.fits'
        cpy_mask['logfile'] = 'ctcubemask_py3.log'
        cpy_mask['chatter'] = 4
        cpy_mask.execute()

        # Check result file
        self._check_result_file('ctcubemask_py3.fits', events=573)

        # Clear ctcubemask tool
        cpy_mask.clear()

        # Check that empty ctcubemask tool holds an empty observation
        self._check_obs(cpy_mask.obs(), nobs=0)

        # Prepare observation container
        obs = self._obs_mixed()
        obs.models(gammalib.GModels(self._model))

        # Set-up ctcubemask from observation container, don't perform any
        # energy selection
        mask = ctools.ctcubemask(obs)
        mask['regfile'] = 'NONE'
        mask['ra']      = 83.63
        mask['dec']     = 22.01
        mask['rad']     = 3.0
        mask['emin']    = 'NONE'
        mask['emax']    = 'NONE'
        mask['outcube'] = 'ctcubemask_py4.xml'
        mask['logfile'] = 'ctcubemask_py4.log'
        mask['chatter'] = 3

        # Execute ctcubemask tool
        mask.logFileOpen()   # Make sure we get a log file
        mask.execute()

        # Check result file
        self._check_result_file('filtered_crab_cntmap.fits', events=5542)

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, events=4921):
        """
        Check result file
        """
        # Load counts cube
        cube = gammalib.GCTAEventCube(filename)

        # Check counts cube
        self.test_value(cube.size(), 800000, 'Check for number of cube bins')
        self.test_value(cube.number(), events, 'Check for number of events')

        # Return
        return

    # Check observation container
    def _check_obs(self, obs, nobs=1):
        """
        Check observation container

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Observation container
        nobs : int, optional
            Expected number of observations
        """
        # Check size of observation container
        self.test_value(obs.size(), nobs, 'Check size of observation container')

        # Return
        return
