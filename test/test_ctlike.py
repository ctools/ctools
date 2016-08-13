#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctlike tool.
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


# ========================== #
# Test class for ctlike tool #
# ========================== #
class Test(test):
    """
    Test class for ctlike tool
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
        self.name('ctlike')

        # Append tests
        self.append(self._test_cmd, 'Test ctlike on command line')
        self.append(self._test_python, 'Test ctlike from Python')

        # Return
        return

    # Test ctlike on command line
    def _test_cmd(self):
        """
        Test ctlike on the command line
        """
        # Set tool name
        ctlike = self._tool('ctlike')

        # Setup ctlike command
        cmd = ctlike+' inobs="'+self._events+'"'+ \
                     ' inmodel="'+self._model+'"'+ \
                     ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                     ' outmodel="ctlike_cmd1.xml"'+ \
                     ' logfile="ctlike_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctlike_cmd1.xml')

        # Setup ctlike command
        cmd = ctlike+' inobs="event_file_that_does_not_exist.fits"'+ \
                     ' inmodel="'+self._model+'"'+ \
                     ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                     ' outmodel="ctlike_cmd2.xml"'+ \
                     ' logfile="ctlike_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctlike --help option
        cmd = ctlike+' --help'

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

    # Test ctlike from Python
    def _test_python(self):
        """
        Test ctlike from Python
        """
        # Allocate ctlike
        like = ctools.ctlike()

        # Check that empty ctlike tool holds an empty observation and
        # optimizer
        self._check_obs(like.obs(), nobs=0, nmodels=0)
        self._check_opt(like.opt())

        # Check that saving saves an empty model definition file
        like['outmodel'] = 'ctlike_py0.xml'
        like['logfile']  = 'ctlike_py0.log'
        like.logFileOpen()
        like.save()
        self._check_result_file('ctlike_py0.xml', nmodels=0)

        # Check that clearing does not lead to an exception or segfault
        like.clear()

        # Now set ctlike parameters
        like['inobs']    = self._events
        like['inmodel']  = self._model
        like['caldb']    = self._caldb
        like['irf']      = self._irf
        like['outmodel'] = 'ctlike_py1.xml'
        like['logfile']  = 'ctlike_py1.log'
        like['chatter']  = 2

        # Run ctlike tool
        like.logFileOpen()   # Make sure we get a log file
        like.run()
        like.save()

        # Check result file
        self._check_result_file('ctlike_py1.xml')

        # Copy ctlike tool
        cpy_like = like.copy()

        # Check observation and optimizer of copy
        self._check_obs(cpy_like.obs())
        self._check_opt(cpy_like.opt())

        # Execute copy of ctlike tool again, now with a higher chatter
        # level than before
        cpy_like['outmodel'] = 'ctlike_py2.xml'
        cpy_like['logfile']  = 'ctlike_py2.log'
        cpy_like['chatter']  = 3
        cpy_like.logFileOpen()  # Needed to get a new log file
        cpy_like.execute()

        # Check result file
        self._check_result_file('ctlike_py2.xml')

        # Now clear copy of ctlike tool
        cpy_like.clear()

        # Check that the cleared copy has also cleared the observations
        self._check_obs(cpy_like.obs(), nobs=0, nmodels=0)

        # Prepare observation container with a single event list
        obs = gammalib.GObservations()
        obs.append(gammalib.GCTAObservation(self._events))

        # Prepare a special model to enforce TS computation
        models = gammalib.GModels(self._model)
        models['Crab'].tscalc(True)
        obs.models(models)

        # Allocate ctlike tool from observation container
        like = ctools.ctlike(obs)
        like['caldb']           = self._caldb
        like['irf']             = self._irf
        like['refit']           = True
        like['fix_spat_for_ts'] = True
        like['outmodel']        = 'ctlike_py3.xml'
        like['logfile']         = 'ctlike_py3.log'
        like['chatter']         = 4

        # Execute ctlike
        like.logFileOpen()  # Needed to get a new log file
        like.execute()

        # Check result file
        self._check_result_file('ctlike_py3.xml')

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
            self.test_value(models['Crab']['Prefactor'].error(), 0.0526982e-16,
                            1.0e-3, 'Check Crab Prefactor error')
            self.test_value(models['Crab']['Index'].value(), -2.43549,
                            1.0e-3, 'Check fitted Crab Index')
            self.test_value(models['Crab']['Index'].error(), 0.0248116,
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
