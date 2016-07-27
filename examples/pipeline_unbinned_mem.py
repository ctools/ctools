#! /usr/bin/env python
# ==========================================================================
# Perform unbinned in-memory analysis of simulated CTA data
#
# Copyright (C) 2015-2016 Juergen Knoedlseder
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
import sys
import gammalib
import ctools
from cscripts import obsutils


# ================== #
# Setup observations #
# ================== #
def setup_observations(pattern='four', ra=83.63, dec=22.01, offset=1.5,
                       emin=0.1, emax=100.0, rad=5.0, duration=180.0,
                       deadc=0.95, caldb='prod2', irf='South_0.5h'):
    """
    Returns an observation container

    Parameters
    ----------
    pattern : str, optional
        Pointing pattern, either 'single' or 'four'
    ra : float, optional
        Right Ascension of pattern centre (deg)
    dec : float, optional
        Declination of pattern centre (deg)
    offset : float, optional
        Offset between observations of pattern (deg)
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    rad : float, optional
        ROI radius used for analysis (deg)
    duration : float, optional
        Duration of one CTA observation (s)
    deadc : float, optional
        Deadtime correction factor
    caldb : str, optional
        Calibration database path
    irf : str, optional
        Instrument response function

    Returns
    -------
    obs : `~gammalib.GObservations`
        Observation container
    """
    # Set list of observations
    obs_def_list = obsutils.set_obs_patterns(pattern,
                                             ra=ra,
                                             dec=dec,
                                             offset=offset)

    # Get observation container
    obs = obsutils.set_obs_list(obs_def_list,
                                duration=duration,
                                emin=emin,
                                emax=emax,
                                rad=rad,
                                caldb=caldb,
                                irf=irf)

    # Return observation container
    return obs


# =========== #
# Setup model #
# =========== #
def setup_model(obs, model='data/crab.xml'):
    """
    Setup model for analysis

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    model : str, optional
        Model definition XML file

    Returns
    -------
    obs : `~gammalib.GObservations`
        Observation container
    """
    # Append model from file to observation container
    obs.models(gammalib.GModels(model))

    # Return observation container
    return obs


# ================================ #
# Simulation and analysis pipeline #
# ================================ #
def run_pipeline(obs, ra=83.63, dec=22.01, rad=3.0,
                 emin=0.1, emax=100.0,
                 tmin=0.0, tmax=0.0,
                 debug=False):
    """
    Simulation and binned analysis pipeline

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    ra : float, optional
        Right Ascension of Region of Interest centre (deg)
    dec : float, optional
        Declination of Region of Interest centre (deg)
    rad : float, optional
        Radius of Region of Interest (deg)
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    tmin : float, optional
        Start time (s)
    tmax : float, optional
        Stop time (s)
    debug : bool, optional
        Debug function
    """
    # Simulate events
    sim = ctools.ctobssim(obs)
    sim['debug'] = debug
    sim.run()

    # Select events
    select = ctools.ctselect(sim.obs())
    select['ra']    = ra
    select['dec']   = dec
    select['rad']   = rad
    select['emin']  = emin
    select['emax']  = emax
    select['tmin']  = tmin
    select['tmax']  = tmax
    select['debug'] = debug
    select.run()

    # Perform maximum likelihood fitting
    like = ctools.ctlike(select.obs())
    like['debug'] = True # Switch this always on for results in console
    like.run()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Dump header
    print('********************************************')
    print('*      CTA unbinned analysis pipeline      *')
    print('********************************************')

    # Get optional arguments
    if len(sys.argv) == 2:
        datadir = sys.argv[1]
    else:
        datadir = 'data'

    # Setup observations
    obs = setup_observations()

    # Setup model
    obs = setup_model(obs, model=datadir+'/crab.xml')

    # Run analysis pipeline
    run_pipeline(obs)
