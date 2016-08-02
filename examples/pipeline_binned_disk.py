#! /usr/bin/env python
# ==========================================================================
# Perform binned analysis of simulated CTA data
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
import cscripts


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
    Simulation and binned analysis pipeline

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    enumbins : int, optional
        Number of energy bins
    nxpix : int, optional
        Number of pixels in X axis
    nypix : int, optional
        Number of pixels in Y axis
    binsz : float, optional
        Pixel size (deg)
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


# ===================================================== #
# Run binned pipeline with intermediate results on disk #
# ===================================================== #
def pipeline_binned_disk():
    """
    Run binned pipeline with intermediate results on disk
    """
    # Set usage string
    usage = 'pipeline_binned_disk.py [-d datadir]'

    # Set default options
    options = [{'option': '-d', 'value': 'data'}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    datadir = options[0]['value']

    # Setup observations
    obs = cscripts.obsutils.set_observations(83.63, 22.01, 5.0, 0.0, 180.0,
                                             0.1, 100.0, 'South_0.5h', 'prod2',
                                             pattern='four', offset=1.5)

    # Setup model
    obs.models(gammalib.GModels(datadir+'/crab.xml'))

    # Run analysis pipeline
    run_pipeline(obs, model=datadir+'/crab.xml', enumbins=10, nxpix=40,
                 nypix=40, binsz=0.1)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Run binned pipeline with intermediate results on disk
    pipeline_binned_disk()
