#! /usr/bin/env python
# ==========================================================================
# Perform binned in-memory analysis of simulated CTA data.
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
                 coordsys='CEL', proj='CAR', debug=False):
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
    debug : bool, optional
        Debug function
    """
    # Simulate events
    sim = ctools.ctobssim(obs)
    sim['debug'] = debug
    sim.run()

    # Bin events by looping over all observations in the container
    obs = gammalib.GObservations()
    obs.models(sim.obs().models())
    for run in sim.obs():

        # Create container with a single observation
        container = gammalib.GObservations()
        container.append(run)

        # Bin events for that observation
        bin = ctools.ctbin(container)
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
        bin.run()

        # Append result to observations
        obs.extend(bin.obs())

    # Perform maximum likelihood fitting
    like = ctools.ctlike(obs)
    like['debug'] = True # Switch this always on for results in console
    like.run()

    # Return
    return


# ============================= #
# Run binned in-memory pipeline #
# ============================= #
def pipeline_binned_mem():
    """
    Run binned in-memory pipeline
    """
    # Set usage string
    usage = 'pipeline_binned_mem.py [-d datadir]'

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
    run_pipeline(obs, enumbins=10, nxpix=40, nypix=40, binsz=0.1)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Run binned in-memory pipeline
    pipeline_binned_mem()
