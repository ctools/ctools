#! /usr/bin/env python
# ==========================================================================
# Perform binned CTA analysis based of simulated CTA data.
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
    Returns an observation container.

    Parameters
    ----------
    pattern : str, optional
        Pointing pattern, either 'single' or 'four'
    ra : float, optional
        Right Ascension of pattern centre [deg]
    dec : float, optional
        Delication of pattern centre [deg]
    offset : float, optional
        Offset between observations of pattern [deg]
    emin : float, optional
        Minimum energy [TeV]
    emax : float, optional
        Maximum energy [TeV]
    rad : float, optional
        ROI radius used for analysis [deg]
    duration : float, optional
        Duration of one CTA observation [seconds]
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
    Setup model for analysis.

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
def run_pipeline(obs, emin=0.1, emax=100.0,
                 enumbins=20, nxpix=200, nypix=200, binsz=0.02,
                 coordsys='CEL', proj='CAR',
                 model='data/crab.xml',
                 caldb='prod2', irf='South_0.5h',
                 debug=False):
    """
    Simulation and binned analysis pipeline.

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    emin : float, optional
        Minimum energy [TeV]
    emax : float, optional
        Maximum energy [TeV]
    enumbins : int, optional
        Number of energy bins
    nxpix : int, optional
        Number of pixels in X axis
    nypix : int, optional
        Number of pixels in Y axis
    binsz : float, optional
        Pixel size [deg/pixel]
    coordsys : str, optional
        Coordinate system
    proj : str, optional
        Coordinate projection
    model : str, optional
        Model definition XML file
    caldb : str, optional
        Calibration database path
    irf : str, optional
        Instrument response function
    debug : bool, optional
        Debug function
    """
    # Simulate events
    sim = ctools.ctobssim(obs)
    sim['debug']     = debug
    sim['outevents'] = 'obs.xml'
    sim.execute()

    # Bin events by looping over all observations in the container
    sim_obs = gammalib.GObservations('obs.xml')
    obs     = gammalib.GObservations()
    for run in sim_obs:

        # Get event filename and set counts cube filename
        eventfile = run.eventfile().url()
        cubefile  = 'cube_'+eventfile

        # Bin events for that observation
        bin = ctools.ctbin()
        bin['inobs']    = eventfile
        bin['outcube']  = cubefile
        bin['ebinalg']  = 'LOG'
        bin['emin']     = emin
        bin['emax']     = emax
        bin['enumbins'] = enumbins
        bin['nxpix']    = nxpix
        bin['nypix']    = nypix
        bin['binsz']    = binsz
        bin['coordsys'] = coordsys
        bin['usepnt']   = True
        bin['proj']     = proj
        bin.execute()

        # Set observation ID
        bin.obs()[0].id(cubefile)
        bin.obs()[0].eventfile(cubefile)

        # Append result to observations
        obs.extend(bin.obs())

    # Save XML file
    xml = gammalib.GXml()
    obs.write(xml)
    xml.save('obs_cube.xml')

    # Perform maximum likelihood fitting
    like = ctools.ctlike()
    like['inobs']    = 'obs_cube.xml'
    like['inmodel']  = model
    like['outmodel'] = 'fit_results.xml'
    like['expcube']  = 'NONE'
    like['psfcube']  = 'NONE'
    like['bkgcube']  = 'NONE'
    like['caldb']    = caldb
    like['irf']      = irf
    like['debug']    = True # Switch this always on for results in console
    like.execute()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Dump header
    print('********************************************')
    print('*       CTA binned analysis pipeline       *')
    print('********************************************')

    # Setup observations
    obs = setup_observations()

    # Setup model
    obs = setup_model(obs)

    # Run analysis pipeline
    run_pipeline(obs)
