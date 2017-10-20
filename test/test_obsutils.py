#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the obsutils module
#
# Copyright (C) 2017 Juergen Knoedlseder
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
from cscripts import obsutils
from testing import test


# =============================== #
# Test class for csviscube script #
# =============================== #
class Test(test):
    """
    Test class for csviscube script
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
        self.name('obsutils')

        # Append tests
        self.append(self._test_sim_unbinned,
                    'Test obsutils.sim() function in unbinned mode')
        self.append(self._test_sim_binned,
                    'Test obsutils.sim() function in binned mode')
        self.append(self._test_sim_stacked,
                    'Test obsutils.sim() function in stacked mode')
        self.append(self._test_sim_stacked_edisp,
                    'Test obsutils.sim() function in stacked mode with energy '
                    'dispersion')
        self.append(self._test_sim_log,
                    'Test obsutils.sim() function with logging switched on')
        self.append(self._test_set_obs,
                    'Test obsutils.set_obs() function')
        self.append(self._test_set_obs_list,
                    'Test obsutils.set_obs_list() function')
        self.append(self._test_set_obs_patterns,
                    'Test obsutils.set_obs_patterns() function')
        self.append(self._test_set_observations,
                    'Test obsutils.set_observations() function')
        self.append(self._test_get_stacked_response,
                    'Test obsutils.get_stacked_response() function')
        self.append(self._test_get_stacked_obs,
                    'Test obsutils.get_stacked_obs() function')

        # Return
        return

    # Setup method for sim() function test
    def _setup_sim(self, two=False):
        """
        Setup method for sim() function test
        """
        # Set-up observation container
        dir = gammalib.GSkyDir()
        dir.radec_deg(83.6331, 22.0145)
        obs = gammalib.GObservations()
        run = obsutils.set_obs(dir, duration=100.0, emin=1.0, emax=10.0)
        run.id('0')
        obs.append(run)
        if two:
            run.id('1')
            obs.append(run)

        # Append model
        obs.models(gammalib.GModels(self._model))

        # Return
        return obs

    # Test sim() function in unbinned mode
    def _test_sim_unbinned(self):
        """
        Test sim() function in unbinned mode
        """
        # Simulate unbinned observations
        res = obsutils.sim(self._setup_sim())

        # Check simulation results
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 2, 'Check number of models')
        self.test_value(res.nobserved(), 40, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'EventList', 'Check event type')
        self.test_value(res[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().emax().TeV(), 10.0,
                        'Check minimum energy')
        self.test_value(res[0].events().number(), 40,
                        'Check number of events in list')

        # Check energy dispersion flag
        self.test_assert(not (res[0].response().use_edisp()),
                         'Check energy dispersion usage')
        res[0].response().apply_edisp(True)
        self.test_assert(res[0].response().use_edisp(),
                         'Check energy dispersion usage')

        # Return
        return

    # Test sim() function in binned mode
    def _test_sim_binned(self):
        """
        Test sim() function in binned mode
        """
        # Simulate binned observations
        res = obsutils.sim(self._setup_sim(), nbins=5)

        # Check simulation results
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 2, 'Check number of models')
        self.test_value(res.nobserved(), 40, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'CountsCube', 'Check event type')
        self.test_value(res[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().emax().TeV(), 10.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().size(), 5,
                        'Check number of energy bins')
        self.test_value(res[0].events().number(), 40,
                        'Check number of events in cube')

        # Check energy dispersion flag
        self.test_assert(not (res[0].response().use_edisp()),
                         'Check energy dispersion usage')
        res[0].response().apply_edisp(True)
        self.test_assert(res[0].response().use_edisp(),
                         'Check energy dispersion usage')

        # Return
        return

    # Test sim() function in stacked mode
    def _test_sim_stacked(self):
        """
        Test sim() function in stacked mode
        """
        # Simulate stacked observations
        res = obsutils.sim(self._setup_sim(two=True), nbins=5)

        # Check simulation results
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 2, 'Check number of models')
        self.test_value(res.nobserved(), 86, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'CountsCube', 'Check event type')
        self.test_value(res[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().emax().TeV(), 10.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().size(), 5,
                        'Check number of energy bins')
        self.test_value(res[0].events().number(), 86,
                        'Check number of events in cube')

        # Check energy dispersion flag
        self.test_assert(not (res[0].response().use_edisp()),
                         'Check energy dispersion usage')
        res[0].response().apply_edisp(True)
        self.test_assert(not (res[0].response().use_edisp()),
                         'Check energy dispersion usage')

        # Return
        return

    # Test sim() function in stacked mode with energy dispersion
    def _test_sim_stacked_edisp(self):
        """
        Test sim() function in stacked mode with energy dispersion
        """
        # Simulate binned observations
        res = obsutils.sim(self._setup_sim(two=True), nbins=5, edisp=True,
                           log=True)

        # Check simulation results
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 2, 'Check number of models')
        self.test_value(res.nobserved(), 82, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'CountsCube', 'Check event type')
        self.test_value(res[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().emax().TeV(), 10.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().size(), 5,
                        'Check number of energy bins')
        self.test_value(res[0].events().number(), 82,
                        'Check number of events in cube')

        # Check energy dispersion flag
        self.test_assert(not (res[0].response().use_edisp()),
                         'Check energy dispersion usage')
        res[0].response().apply_edisp(True)
        self.test_assert(res[0].response().use_edisp(),
                         'Check energy dispersion usage')

        # Return
        return

    # Test sim() function with logging switched on
    def _test_sim_log(self):
        """
        Test sim() function with logging switched on
        """
        # Simulate unbinned observations
        res = obsutils.sim(self._setup_sim(), log=True)

        # Check simulation results
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 2, 'Check number of models')
        self.test_value(res.nobserved(), 40, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'EventList', 'Check event type')
        self.test_value(res[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().emax().TeV(), 10.0,
                        'Check minimum energy')

        # Check energy dispersion flag
        self.test_assert(not (res[0].response().use_edisp()),
                         'Check energy dispersion usage')
        res[0].response().apply_edisp(True)
        self.test_assert(res[0].response().use_edisp(),
                         'Check energy dispersion usage')

        # Return
        return

    # Test set_obs() function
    def _test_set_obs(self):
        """
        Test set_obs() function
        """
        # TODO: implement
        # Return
        return

    # Test set_obs_list() function
    def _test_set_obs_list(self):
        """
        Test set_obs_list() function
        """
        # TODO: implement
        # Return
        return

    # Test set_obs_patterns() function
    def _test_set_obs_patterns(self):
        """
        Test set_obs_pattern() function
        """
        # Get single observation pattern
        res = obsutils.set_obs_patterns('single')
        self.test_value(len(res), 1, 'Check number of observations for '
                        '"single" pattern')

        # Get four observation pattern
        res = obsutils.set_obs_patterns('four')
        self.test_value(len(res), 4, 'Check number of observations for "four"'
                        ' pattern')

        # Get invalud observations pattern
        self.test_try('Check invalid pointing pattern')
        try:
            res = obsutils.set_obs_patterns('five')
            self.test_try_failure('Exception not thrown')
        except RuntimeError:
            self.test_try_success()

        # Return
        return

    # Test set_observations() function
    def _test_set_observations(self):
        """
        Test set_observations() function
        """
        # TODO: implement
        # Return
        return

    # Test get_stacked_response() function
    def _test_get_stacked_response(self):
        """
        Test get_stacked_response() function
        """
        # TODO: implement
        # Return
        return

    # Test get_stacked_obs() function
    def _test_get_stacked_obs(self):
        """
        Test get_stacked_obs() function
        """
        # TODO: implement
        # Return
        return

