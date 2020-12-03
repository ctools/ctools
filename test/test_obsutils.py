#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the obsutils module
#
# Copyright (C) 2017-2020 Juergen Knoedlseder
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
from cscripts import obsutils
from testing import test


# ============================== #
# Test class for obsutils module #
# ============================== #
class Test(test):
    """
    Test class for obsutils module
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set test data
        self._model = self._datadir + '/model_crab_radialacceptance.xml'

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
        self.append(self._test_sim_onoff,
                    'Test obsutils.sim() function in On/Off mode')
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
        self.append(self._test_get_onoff_obs,
                    'Test obsutils.get_onoff_obs() function')

        # Return
        return

    # Setup method for sim() function test
    def _setup_sim(self, two=False):
        """
        Setup method for sim() function test
        """
        # Set-up observation container
        pnt = gammalib.GSkyDir()
        pnt.radec_deg(83.6331, 22.0145)
        obs = gammalib.GObservations()
        run = obsutils.set_obs(pnt, duration=20.0, emin=1.0, emax=10.0)
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
        self.test_value(res.nobserved(), 4, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'EventList', 'Check event type')
        self.test_value(res[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy')
        self.test_value(res[0].events().number(), 4,
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
        self.test_value(res.nobserved(), 4, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'CountsCube', 'Check event type')
        self.test_value(res[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy')
        self.test_value(res[0].events().ebounds().size(), 5,
                        'Check number of energy bins')
        self.test_value(res[0].events().number(), 4,
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
        self.test_value(res.nobserved(), 15, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'CountsCube', 'Check event type')
        self.test_value(res[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy')
        self.test_value(res[0].events().ebounds().size(), 5,
                        'Check number of energy bins')
        self.test_value(res[0].events().number(), 15,
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
        self.test_value(res.nobserved(), 20, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'CountsCube', 'Check event type')
        self.test_value(res[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy')
        self.test_value(res[0].events().ebounds().size(), 5,
                        'Check number of energy bins')
        self.test_value(res[0].events().number(), 20,
                        'Check number of events in cube')

        # Check energy dispersion flag
        self.test_assert(not (res[0].response().use_edisp()),
                         'Check energy dispersion usage')
        res[0].response().apply_edisp(True)
        self.test_assert(res[0].response().use_edisp(),
                         'Check energy dispersion usage')

        # Return
        return

    # Test sim() function in On/Off mode
    def _test_sim_onoff(self):
        """
        Test sim() function in On/Off mode
        """
        # Set-up observation container
        pnt = gammalib.GSkyDir()
        pnt.radec_deg(83.63, 22.51)
        obs = gammalib.GObservations()
        run = obsutils.set_obs(pnt, duration=100.0, emin=1.0, emax=10.0, obsid='0')
        obs.append(run)
        pnt.radec_deg(83.63, 21.51)
        run = obsutils.set_obs(pnt, duration=100.0, emin=1.0, emax=10.0, obsid='1')
        obs.append(run)
        obs.models(gammalib.GModels(self._model))

        # Simulate On/Off observations
        res = obsutils.sim(obs, onsrc='Crab', nbins=5)

        # Check simulation results
        self.test_value(res.size(), 2, 'Check number of observations')
        self.test_value(res.models().size(), 2, 'Check number of models')
        self.test_value(res.nobserved(), 48, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')

        # Check results of first observation
        self.test_value(res[0].on_spec().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy of On spectrum')
        self.test_value(res[0].on_spec().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy of On spectrum')
        self.test_value(res[0].on_spec().ebounds().size(), 5,
                        'Check number of energy bins of On spectrum')
        self.test_value(res[0].on_spec().counts(), 26,
                        'Check number of events in of On spectrum')
        self.test_value(res[0].off_spec().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy of Off spectrum')
        self.test_value(res[0].off_spec().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy of Off spectrum')
        self.test_value(res[0].off_spec().ebounds().size(), 5,
                        'Check number of energy bins of Off spectrum')
        self.test_value(res[0].off_spec().counts(), 1,
                        'Check number of events in of Off spectrum')
        self.test_value(res[0].arf().ebounds().emin().TeV(), 0.5,
                        'Check minimum energy of ARF')
        self.test_value(res[0].arf().ebounds().emax().TeV(), 12.0,
                        'Check maximum energy of ARF')
        self.test_value(res[0].arf().ebounds().size(), 41,
                        'Check number of energy bins of ARF')
        self.test_value(res[0].rmf().etrue().emin().TeV(), 0.5,
                        'Check minimum true energy of RMF')
        self.test_value(res[0].rmf().etrue().emax().TeV(), 12.0,
                        'Check maximum true energy of RMF')
        self.test_value(res[0].rmf().etrue().size(), 41,
                        'Check number of true energy bins of RMF')
        self.test_value(res[0].rmf().emeasured().emin().TeV(), 1.0,
                        'Check minimum reconstructed energy of RMF')
        self.test_value(res[0].rmf().emeasured().emax().TeV(), 10.0,
                        'Check maximum reconstructed energy of RMF')
        self.test_value(res[0].rmf().emeasured().size(), 5,
                        'Check number of reconstructed energy bins of RMF')

        # Check results of second observation
        self.test_value(res[1].on_spec().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy of On spectrum')
        self.test_value(res[1].on_spec().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy of On spectrum')
        self.test_value(res[1].on_spec().ebounds().size(), 5,
                        'Check number of energy bins of On spectrum')
        self.test_value(res[1].on_spec().counts(), 22,
                        'Check number of events in of On spectrum')
        self.test_value(res[1].off_spec().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy of Off spectrum')
        self.test_value(res[1].off_spec().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy of Off spectrum')
        self.test_value(res[1].off_spec().ebounds().size(), 5,
                        'Check number of energy bins of Off spectrum')
        self.test_value(res[1].off_spec().counts(), 1,
                        'Check number of events in of Off spectrum')
        self.test_value(res[1].arf().ebounds().emin().TeV(), 0.5,
                        'Check minimum energy of ARF')
        self.test_value(res[1].arf().ebounds().emax().TeV(), 12.0,
                        'Check maximum energy of ARF')
        self.test_value(res[1].arf().ebounds().size(), 41,
                        'Check number of energy bins of ARF')
        self.test_value(res[1].rmf().etrue().emin().TeV(), 0.5,
                        'Check minimum true energy of RMF')
        self.test_value(res[1].rmf().etrue().emax().TeV(), 12.0,
                        'Check maximum true energy of RMF')
        self.test_value(res[1].rmf().etrue().size(), 41,
                        'Check number of true energy bins of RMF')
        self.test_value(res[1].rmf().emeasured().emin().TeV(), 1.0,
                        'Check minimum reconstructed energy of RMF')
        self.test_value(res[1].rmf().emeasured().emax().TeV(), 10.0,
                        'Check maximum reconstructed energy of RMF')
        self.test_value(res[1].rmf().emeasured().size(), 5,
                        'Check number of reconstructed energy bins of RMF')

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
        self.test_value(res.nobserved(), 4, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'EventList', 'Check event type')
        self.test_value(res[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy')

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
        # Setup pointing direction
        pnt = gammalib.GSkyDir()
        pnt.radec_deg(83.63, 22.51)

        # Setup one CTA observation
        res = obsutils.set_obs(pnt, emin=1.0, emax=10.0)

        # Check result
        self.test_value(res.eventtype(), 'EventList', 'Check event type')
        self.test_value(res.events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res.events().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy')
        self.test_value(res.pointing().dir().ra_deg(), 83.63,
                        'Check pointing Right Ascension')
        self.test_value(res.pointing().dir().dec_deg(), 22.51,
                        'Check pointing declination')

        # Setup one CTA observation for local caldb directory
        res = obsutils.set_obs(pnt, caldb=self._datadir, irf='irf_file.fits',
                               emin=1.0, emax=10.0)

        # Check result
        self.test_value(res.eventtype(), 'EventList', 'Check event type')
        self.test_value(res.events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res.events().ebounds().emax().TeV(), 10.0,
                        'Check maximum energy')
        self.test_value(res.pointing().dir().ra_deg(), 83.63,
                        'Check pointing Right Ascension')
        self.test_value(res.pointing().dir().dec_deg(), 22.51,
                        'Check pointing declination')

        # Return
        return

    # Test set_obs_list() function
    def _test_set_obs_list(self):
        """
        Test set_obs_list() function
        """
        # Setup observation definition list
        obsdeflist = [{'ra': 83.63, 'dec': 22.51}]

        # Setup observation container
        res = obsutils.set_obs_list(obsdeflist)

        # Check result
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 0, 'Check number of models')
        self.test_value(res.nobserved(), 0, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'EventList', 'Check event type')
        self.test_value(res[0].pointing().dir().ra_deg(), 83.63,
                        'Check pointing Right Ascension')
        self.test_value(res[0].pointing().dir().dec_deg(), 22.51,
                        'Check pointing declination')

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

        # Get invalid observations pattern
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
        # Set observation
        res = obsutils.set_observations(83.63, 22.51, 5.0,
                                        0.0, 1800.0,
                                        1.0, 100.0,
                                        'South_50h', 'prod2',
                                        pattern='single')

        # Check result
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 0, 'Check number of models')
        self.test_value(res.nobserved(), 0, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_value(res[0].eventtype(), 'EventList', 'Check event type')
        self.test_value(res[0].events().ebounds().emin().TeV(), 1.0,
                        'Check minimum energy')
        self.test_value(res[0].events().ebounds().emax().TeV(), 100.0,
                        'Check maximum energy')
        self.test_value(res[0].pointing().dir().ra_deg(), 83.63,
                        'Check pointing Right Ascension')
        self.test_value(res[0].pointing().dir().dec_deg(), 22.51,
                        'Check pointing declination')

        # Return
        return

    # Test get_stacked_response() function
    def _test_get_stacked_response(self):
        """
        Test get_stacked_response() function
        """
        # Set-up observation container
        obs = self._setup_sim()

        # Set-up counts cube
        map     = gammalib.GSkyMap('CAR','CEL',83.6331,22.0145,0.1,0.1,10,10,5)
        emin    = gammalib.GEnergy(0.1, 'TeV')
        emax    = gammalib.GEnergy(100.0, 'TeV')
        ebds    = gammalib.GEbounds(5, emin, emax)
        tmin    = gammalib.GTime(0.0)
        tmax    = gammalib.GTime(1000.0)
        gti     = gammalib.GGti(tmin, tmax)
        cntcube = gammalib.GCTAEventCube(map, ebds, gti)

        # Get stacked response
        res = obsutils.get_stacked_response(obs, cntcube, edisp=False)

        # Check result
        self.test_value(res['expcube'].cube().npix(), 100,
                        'Check number of exposure cube pixels')
        self.test_value(res['expcube'].cube().nmaps(), 91,
                        'Check number of exposure cube maps')
        self.test_value(res['psfcube'].cube().npix(), 4,
                        'Check number of PSF cube pixels')
        self.test_value(res['psfcube'].cube().nmaps(), 18200,
                        'Check number of PSF cube maps')
        self.test_value(res['bkgcube'].cube().npix(), 100,
                        'Check number of background cube pixels')
        self.test_value(res['bkgcube'].cube().nmaps(), 5,
                        'Check number of background cube maps')

        # Get stacked response with energy dispersion
        res = obsutils.get_stacked_response(obs, cntcube, edisp=True)

        # Check result
        self.test_value(res['expcube'].cube().npix(), 100,
                        'Check number of exposure cube pixels')
        self.test_value(res['expcube'].cube().nmaps(), 105,
                        'Check number of exposure cube maps')
        self.test_value(res['psfcube'].cube().npix(), 4,
                        'Check number of PSF cube pixels')
        self.test_value(res['psfcube'].cube().nmaps(), 21000,
                        'Check number of PSF cube maps')
        self.test_value(res['bkgcube'].cube().npix(), 100,
                        'Check number of background cube pixels')
        self.test_value(res['bkgcube'].cube().nmaps(), 5,
                        'Check number of background cube maps')
        self.test_value(res['edispcube'].cube().npix(), 4,
                        'Check number of energy dispersion cube pixels')
        self.test_value(res['edispcube'].cube().nmaps(), 10500,
                        'Check number of energy dispersion cube maps')

        # Set-up counts cube with large number of energy bins
        ebds    = gammalib.GEbounds(100, emin, emax)
        cntcube = gammalib.GCTAEventCube(map, ebds, gti)

        # Get stacked response with xref/yref set and large number of energy bins
        res = obsutils.get_stacked_response(obs, cntcube, edisp=False)

        # Check result
        self.test_value(res['expcube'].cube().npix(), 100,
                        'Check number of exposure cube pixels')
        self.test_value(res['expcube'].cube().nmaps(), 91,
                        'Check number of exposure cube maps')
        self.test_value(res['psfcube'].cube().npix(), 4,
                        'Check number of PSF cube pixels')
        self.test_value(res['psfcube'].cube().nmaps(), 18200,
                        'Check number of PSF cube maps')
        self.test_value(res['bkgcube'].cube().npix(), 100,
                        'Check number of background cube pixels')
        self.test_value(res['bkgcube'].cube().nmaps(), 100,
                        'Check number of background cube maps')

        # Return
        return

    # Test get_stacked_obs() function
    def _test_get_stacked_obs(self):
        """
        Test get_stacked_obs() function
        """
        # Set-up unbinned cslightcrv (is not run here!!!)
        lcrv = cscripts.cslightcrv()
        lcrv['inobs']    = self._events
        lcrv['inmodel']  = self._model
        lcrv['srcname']  = 'Crab'
        lcrv['caldb']    = self._caldb
        lcrv['irf']      = self._irf
        lcrv['tbinalg']  = 'LIN'
        lcrv['tmin']     = '2020-01-01T00:00:00'
        lcrv['tmax']     = '2020-01-01T00:05:00'
        lcrv['tbins']    = 2
        lcrv['method']   = '3D'
        lcrv['emin']     = 1.0
        lcrv['emax']     = 100.0
        lcrv['enumbins'] = 3
        lcrv['coordsys'] = 'CEL'
        lcrv['proj']     = 'TAN'
        lcrv['xref']     = 83.63
        lcrv['yref']     = 22.01
        lcrv['binsz']    = 0.1
        lcrv['nxpix']    = 10
        lcrv['nypix']    = 10
        lcrv['outfile']  = 'obsutils_cslightcrv_py1.fits'
        lcrv['logfile']  = 'obsutils_cslightcrv_py1.log'
        lcrv._get_parameters()

        # Get stacked observation container
        res = obsutils.get_stacked_obs(lcrv, lcrv.obs())

        # Apply energy dispersion
        res[0].response().apply_edisp(True)

        # Check result
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 2, 'Check number of models')
        self.test_value(res.nobserved(), 96, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_assert(not res[0].response().use_edisp(),
                         'Check that no energy dispersion is used')

        # Set energy dispersion
        lcrv['edisp'] = True

        # Get stacked observation container
        res = obsutils.get_stacked_obs(lcrv, lcrv.obs())

        # Apply energy dispersion
        res[0].response().apply_edisp(True)

        # Check result
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 2, 'Check number of models')
        self.test_value(res.nobserved(), 96, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_assert(res[0].response().use_edisp(),
                         'Check that energy dispersion is used')

        # Increase chatter level
        lcrv['edisp']   = False
        lcrv['chatter'] = 4

        # Get stacked observation container
        res = obsutils.get_stacked_obs(lcrv, lcrv.obs())

        # Apply energy dispersion
        res[0].response().apply_edisp(True)

        # Check result
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 2, 'Check number of models')
        self.test_value(res.nobserved(), 96, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')
        self.test_assert(not res[0].response().use_edisp(),
                         'Check that no energy dispersion is used')

        # Return
        return

    # Test get_onoff_obs() function
    def _test_get_onoff_obs(self):
        """
        Test get_onoff_obs() function
        """
        # Set-up unbinned cslightcrv (is not run here!!!)
        lcrv = cscripts.cslightcrv()
        lcrv['inobs']    = self._events
        lcrv['inmodel']  = self._model
        lcrv['srcname']  = 'Crab'
        lcrv['caldb']    = self._caldb
        lcrv['irf']      = self._irf
        lcrv['tbinalg']  = 'LIN'
        lcrv['tmin']     = '2020-01-01T00:00:00'
        lcrv['tmax']     = '2020-01-01T00:05:00'
        lcrv['tbins']    = 2
        lcrv['method']   = 'ONOFF'
        lcrv['srcshape'] = 'CIRCLE'
        lcrv['emin']     = 1.0
        lcrv['emax']     = 100.0
        lcrv['enumbins'] = 2
        lcrv['coordsys'] = 'CEL'
        lcrv['xref']     = 83.63
        lcrv['yref']     = 22.01
        lcrv['rad']      = 0.2
        lcrv['outfile']  = 'obsutils_cslightcrv_py2.fits'
        lcrv['logfile']  = 'obsutils_cslightcrv_py2.log'
        lcrv._get_parameters()

        # Get on/off observation container
        res = obsutils.get_onoff_obs(lcrv, lcrv.obs())

        # Check result
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 2, 'Check number of models')
        self.test_value(res.nobserved(), 91, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')

        # Use galactic coordinates
        lcrv['coordsys'] = 'GAL'
        lcrv['xref']     = 184.5597
        lcrv['yref']     =  -5.7892
        lcrv['chatter']  = 4

        # Get on/off observation container
        res = obsutils.get_onoff_obs(lcrv, lcrv.obs())

        # Check result
        self.test_value(res.size(), 1, 'Check number of observations')
        self.test_value(res.models().size(), 2, 'Check number of models')
        self.test_value(res.nobserved(), 91, 'Check number of observed events')
        self.test_value(res.npred(), 0.0, 'Check number of predicted events')

        # Return
        return
