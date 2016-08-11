#! /usr/bin/env python
# ==========================================================================
# This scripts performs tests processing pipelines
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


# ============================ #
# Test class for ctobssim tool #
# ============================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for pipelines
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Set test data directory
        self._datadir = os.environ['TEST_DATA']

        # Set some standard test data
        self._model   = self._datadir + '/crab.xml'
        self._caldb   = 'prod2'
        self._irf     = 'South_0.5h'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('pipelines')

        # Append tests
        self.append(self.test_unbinned_fits,
                    'Test unbinned pipeline with FITS file saving')
        self.append(self.test_unbinned_mem,
                    'Test unbinned in-memory pipeline')

        # Return
        return

    # Test unbinned pipeline with FITS file saving
    def test_unbinned_fits(self):
        """
        Test unbinned pipeline with FITS file saving
        """
        # Set script parameters
        events_name          = 'events.fits'
        selected_events_name = 'selected_events.fits'
        result_name          = 'results.xml'
        ra                   =   83.63
        dec                  =   22.01
        rad_sim              =   10.0
        rad_select           =    3.0
        tstart               =    0.0
        tstop                = 1800.0
        emin                 =    0.1
        emax                 =  100.0

        # Simulate events
        sim = ctools.ctobssim()
        sim['inmodel']   = self._model
        sim['outevents'] = events_name
        sim['caldb']     = self._caldb
        sim['irf']       = self._irf
        sim['ra']        = ra
        sim['dec']       = dec
        sim['rad']       = rad_sim
        sim['tmin']      = tstart
        sim['tmax']      = tstop
        sim['emin']      = emin
        sim['emax']      = emax
        sim.execute()

        # Select events
        select = ctools.ctselect()
        select['inobs']  = events_name
        select['outobs'] = selected_events_name
        select['ra']     = ra
        select['dec']    = dec
        select['rad']    = rad_select
        select['tmin']   = tstart
        select['tmax']   = tstop
        select['emin']   = emin
        select['emax']   = emax
        select.execute()

        # Perform maximum likelihood fitting
        like = ctools.ctlike()
        like['inobs']    = selected_events_name
        like['inmodel']  = self._model
        like['outmodel'] = result_name
        like['caldb']    = self._caldb
        like['irf']      = self._irf
        like.execute()

        # Return
        return

    # Test unbinned in-memory pipeline
    def test_unbinned_mem(self):
        """
        Test unbinned in-memory pipeline
        """
        # Set script parameters
        ra         =   83.63
        dec        =   22.01
        rad_sim    =   10.0
        rad_select =    3.0
        tstart     =    0.0
        tstop      = 1800.0
        emin       =    0.1
        emax       =  100.0

        # Simulate events
        sim = ctools.ctobssim()
        sim['inmodel'] = self._model
        sim['caldb']   = self._caldb
        sim['irf']     = self._irf
        sim['ra']      = ra
        sim['dec']     = dec
        sim['rad']     = rad_sim
        sim['tmin']    = tstart
        sim['tmax']    = tstop
        sim['emin']    = emin
        sim['emax']    = emax
        sim.run()

        # Select events
        select = ctools.ctselect(sim.obs())
        select['ra']   = ra
        select['dec']  = dec
        select['rad']  = rad_select
        select['tmin'] = tstart
        select['tmax'] = tstop
        select['emin'] = emin
        select['emax'] = emax
        select.run()

        # Perform maximum likelihood fitting
        like = ctools.ctlike(select.obs())
        like.run()

        # Return
        return
