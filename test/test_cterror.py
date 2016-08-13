#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cterror tool.
#
# Copyright (C) 2015-2016 Florent Forest
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


# =========================== #
# Test class for cterror tool #
# =========================== #
class Test(test):
    """
    Test class for cterror tool
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
        Set all test functions
        """
        # Set test name
        self.name('cterror')

        # Append tests
        self.append(self._test_cmd, 'Test cterror on command line')
        self.append(self._test_python, 'Test cterror from Python')

        # Return
        return

    # Test cterror on command line
    def _test_cmd(self):
        """
        Test cterror on the command line
        """
        # Set tool name
        cterror = self._tool('cterror')

        # Setup cterror command
        cmd = cterror+' inobs="'+self._events+'"'+ \
                      ' inmodel="'+self._model+'" srcname="Crab"'+ \
                      ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                      ' outmodel="cterror_cmd1.xml"'+ \
                      ' logfile="cterror_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('cterror_cmd1.xml')

        # Setup cterror command
        cmd = cterror+' inobs="event_file_that_does_not_exist.fits"'+ \
                      ' inmodel="'+self._model+'" srcname="Crab"'+ \
                      ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                      ' outmodel="cterror_cmd2.xml"'+ \
                      ' logfile="cterror_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup cterror --help option
        cmd = cterror+' --help'

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

    # Test cterror from Python
    def _test_python(self):
        """
        Test cterror from Python
        """
        # Allocate cterror
        error = ctools.cterror()

        # Check that empty cterror tool holds an empty observation and
        # optimizer
        self._check_obs(error.obs(), nobs=0, nmodels=0)
        self._check_opt(error.opt())

        # Check that saving saves an empty model definition file
        error['outmodel'] = 'cterror_py0.xml'
        error['logfile']  = 'cterror_py0.log'
        error.logFileOpen()
        error.save()
        self._check_result_file('cterror_py0.xml', nmodels=0)

        # Check that clearing does not lead to an exception or segfault
        error.clear()

        # Now set cterror parameters
        error['inobs']    = self._events
        error['inmodel']  = self._model
        error['srcname']  = 'Crab'
        error['caldb']    = self._caldb
        error['irf']      = self._irf
        error['outmodel'] = 'cterror_py1.xml'
        error['logfile']  = 'cterror_py1.log'
        error['chatter']  = 2

        # Run cterror tool
        error.logFileOpen()   # Make sure we get a log file
        error.run()
        error.save()

        # Check result file
        self._check_result_file('cterror_py1.xml')

        # Copy cterror tool
        cpy_error = error.copy()

        # Check observation and optimizer of copy
        self._check_obs(cpy_error.obs())
        self._check_opt(cpy_error.opt())

        # Execute copy of cterror tool again, now with a higher chatter
        # level than before
        cpy_error['outmodel'] = 'cterror_py2.xml'
        cpy_error['logfile']  = 'cterror_py2.log'
        cpy_error['chatter']  = 3
        cpy_error.logFileOpen()  # Needed to get a new log file
        cpy_error.execute()

        # Check result file
        self._check_result_file('cterror_py2.xml')

        # Now clear copy of cterror tool
        cpy_error.clear()

        # Check that the cleared copy has also cleared the observations
        self._check_obs(cpy_error.obs(), nobs=0, nmodels=0)

        # Prepare observation container with a single event list
        obs = gammalib.GObservations()
        obs.append(gammalib.GCTAObservation(self._events))
        obs.models(gammalib.GModels(self._model))

        # Setup cterror tool from observation container
        error = ctools.cterror(obs)
        error['srcname']  = 'Crab'
        error['caldb']    = self._caldb
        error['irf']      = self._irf
        error['outmodel'] = 'cterror_py3.xml'
        error['logfile']  = 'cterror_py3.log'
        error['chatter']  = 4

        # Execute cterror tool
        error.logFileOpen()   # Make sure we get a log file
        error.execute()

        # Check result file
        self._check_result_file('cterror_py3.xml')

        # And now a run with not enough iterations
        error['max_iter'] = 1
        error['outmodel'] = 'cterror_py4.xml'
        error['logfile']  = 'cterror_py4.log'
        
        # Execute cterror tool and catch the exception
        error.logFileOpen()   # Make sure we get a log file
        self.test_try('Test cterror with not enough iterations')
        try:
            error.execute()
            self.test_try_failure('Exception not thrown')
        except ValueError:
            self.test_try_success()

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nmodels=2):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Model definition XML file
        nmodels : int, optional
            Expected number of models
        """
        # Open models from result file
        models = gammalib.GModels(filename)

        # Check models
        self._check_models(models, nmodels=nmodels)

        # Return
        return

    # Check observation and models
    def _check_obs(self, obs, nobs=1, nmodels=2):
        """
        Check observation and models

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Models
        nobs : int, optional
            Expected number of observations
        nmodels : int, optional
            Expected number of models
        """
        # Check number of observations
        self.test_value(obs.size(), nobs, 'Check number of observations')

        # Check models
        self._check_models(obs.models(), nmodels=nmodels)

        # Return
        return

    # Check models
    def _check_models(self, models, nmodels=2):
        """
        Check observation and models

        Parameters
        ----------
        models : `~gammalib.GModels`
            Models
        nmodels : int, optional
            Expected number of models
        """
        # Check number of models
        self.test_value(models.size(), nmodels, 'Check number of models')

        # If there are models in the container then check them
        if nmodels > 0:
            self.test_value(models['Crab']['Prefactor'].value(), 1.58907e-16,
                            1.0e-3, 'Check fitted Crab Prefactor')
            self.test_value(models['Crab']['Prefactor'].error(), 0.0529105e-16,
                            1.0e-3, 'Check Crab Prefactor error')
            self.test_value(models['Crab']['Index'].value(), -2.43549,
                            1.0e-3, 'Check fitted Crab Index')
            self.test_value(models['Crab']['Index'].error(), 0.027804,
                            1.0e-3, 'Check Crab Index error')
            self.test_value(models['Background']['Prefactor'].value(), 61.6919e-6,
                            1.0e-3, 'Check fitted background Prefactor')
            self.test_value(models['Background']['Prefactor'].error(), 1.49438e-6,
                            1.0e-3, 'Check background Prefactor error')
            self.test_value(models['Background']['Index'].value(), -2.20535,
                            1.0e-3, 'Check fitted background Index')
            self.test_value(models['Background']['Index'].error(), 0.0113269,
                            1.0e-3, 'Check background Index error')
            self.test_value(models['Background']['Sigma'].value(), 3.04252,
                            1.0e-3, 'Check fitted background Sigma')
            self.test_value(models['Background']['Sigma'].error(), 0.0307008,
                            1.0e-3, 'Check background Sigma error')

        # Return
        return

    # Check optimizer
    def _check_opt(self, opt):
        """
        Check optimizer

        Parameters
        ----------
        opt : `~gammalib.GOptimizerLM`
            Optimizer
        """
        # Check optimizer
        self.test_value(opt.status(), 0, 'Check optimizer status')

        # Return
        return
