#! /usr/bin/env python
# ==========================================================================
# This script illustrates how to perform a stacked CTA analysis based on
# simulated CTA data. You may use and adapt this script to implement your
# own pipeline.
#
# Usage:
# ./pipeline_stacked_disk.py
#
# ==========================================================================
#
# Copyright (C) 2015 Juergen Knoedlseder
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
def setup_observations(pattern="four", ra=83.63, dec=22.01, offset=1.5,
                       emin=0.1, emax=100.0, rad=5.0, duration=180.0,
                       deadc=0.95,
                       caldb="prod2", irf="South_50h"):
    """
    Returns an observation container.

    Keywords:
     pattern   - Pointing pattern, either "single" or "four"
     ra        - RA of pattern centre [deg] (default: 83.6331)
     dec       - DEC of pattern centre [deg] (default: 22.0145)
     offset    - Offset between observations of pattern [deg] (default: 1.5)
     emin      - Minimum energy [TeV] (default: 0.1)
     emax      - Maximum energy [TeV] (default: 100.0)
     rad       - ROI radius used for analysis [deg] (default: 5.0)
     duration  - Duration of one CTA observation [seconds] (default: 1800.0)
     deadc     - Deadtime correction factor (default: 0.95)
     caldb     - Calibration database path (default: "dummy")
     irf       - Instrument response function (default: cta_dummy_irf)
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
def setup_model(obs, model="${CTOOLS}/share/models/crab.xml"):
    """
    Setup model for analysis.

    Keywords:
     model - Model Xml file
    """
    # Append model from file to observation container
    obs.models(gammalib.GModels(model))

    # Return observation container
    return obs


# ================================ #
# Simulation and analysis pipeline #
# ================================ #
def run_pipeline(obs, ra=83.63, dec=22.01, emin=0.1, emax=100.0,
                 enumbins=20, nxpix=200, nypix=200, binsz=0.02,
                 coordsys="CEL", proj="CAR",
                 model="${CTOOLS}/share/models/crab.xml",
                 caldb="prod2", irf="South_50h",
                 debug=False):
    """
    Simulation and stacked analysis pipeline.

    Keywords:
     ra       - RA of cube centre [deg] (default: 83.6331)
     dec      - DEC of cube centre [deg] (default: 22.0145)
     emin     - Minimum energy of cube [TeV] (default: 0.1)
     emax     - Maximum energy of cube [TeV] (default: 100.0)
     enumbins - Number of energy bins in cube (default: 20)
     nxpix    - Number of RA pixels in cube (default: 200)
     nypix    - Number of DEC pixels in cube (default: 200)
     binsz    - Spatial cube bin size [deg] (default: 0.02)
     coordsys - Cube coordinate system (CEL or GAL)
     proj     - Cube World Coordinate System (WCS) projection
     debug    - Enable debugging (default: False)
    """
    # Simulate events
    sim = ctools.ctobssim(obs)
    sim["debug"]     = debug
    sim["outevents"] = "obs.xml"
    sim.execute()

    # Bin events into counts map
    bin = ctools.ctbin()
    bin["inobs"]    = "obs.xml"
    bin["outcube"]  = "cntcube.fits"
    bin["ebinalg"]  = "LOG"
    bin["emin"]     = emin
    bin["emax"]     = emax
    bin["enumbins"] = enumbins
    bin["nxpix"]    = nxpix
    bin["nypix"]    = nypix
    bin["binsz"]    = binsz
    bin["coordsys"] = coordsys
    bin["proj"]     = proj
    bin["xref"]     = ra
    bin["yref"]     = dec
    bin["debug"]    = debug
    bin.execute()

    # Create exposure cube
    expcube = ctools.ctexpcube()
    expcube["inobs"]    = "obs.xml"
    expcube["incube"]   = "cntcube.fits"
    expcube["outcube"]  = "expcube.fits"
    expcube["caldb"]    = caldb
    expcube["irf"]      = irf
    expcube["ebinalg"]  = "LOG"
    expcube["emin"]     = emin
    expcube["emax"]     = emax
    expcube["enumbins"] = enumbins
    expcube["nxpix"]    = nxpix
    expcube["nypix"]    = nypix
    expcube["binsz"]    = binsz
    expcube["coordsys"] = coordsys
    expcube["proj"]     = proj
    expcube["xref"]     = ra
    expcube["yref"]     = dec
    expcube["debug"]    = debug
    expcube.execute()

    # Create PSF cube
    psfcube = ctools.ctpsfcube()
    psfcube["inobs"]    = "obs.xml"
    psfcube["incube"]   = "NONE"
    psfcube["outcube"]  = "psfcube.fits"
    psfcube["caldb"]    = caldb
    psfcube["irf"]      = irf
    psfcube["ebinalg"]  = "LOG"
    psfcube["emin"]     = emin
    psfcube["emax"]     = emax
    psfcube["enumbins"] = enumbins
    psfcube["nxpix"]    = 10
    psfcube["nypix"]    = 10
    psfcube["binsz"]    = 1.0
    psfcube["coordsys"] = coordsys
    psfcube["proj"]     = proj
    psfcube["xref"]     = ra
    psfcube["yref"]     = dec
    psfcube["debug"]    = debug
    psfcube.execute()

    # Create background cube
    bkgcube = ctools.ctbkgcube()
    bkgcube["inobs"]    = "obs.xml"
    bkgcube["inmodel"]  = model
    bkgcube["incube"]   = "cntcube.fits"
    bkgcube["outcube"]  = "bkgcube.fits"
    bkgcube["outmodel"] = "model_bkg.xml"
    bkgcube["caldb"]    = caldb
    bkgcube["irf"]      = irf
    bkgcube["ebinalg"]  = "LOG"
    bkgcube["emin"]     = emin
    bkgcube["emax"]     = emax
    bkgcube["enumbins"] = enumbins
    bkgcube["nxpix"]    = 10
    bkgcube["nypix"]    = 10
    bkgcube["binsz"]    = 1.0
    bkgcube["coordsys"] = coordsys
    bkgcube["proj"]     = proj
    bkgcube["xref"]     = ra
    bkgcube["yref"]     = dec
    bkgcube["debug"]    = debug
    bkgcube.execute()

    # Perform maximum likelihood fitting
    like = ctools.ctlike()
    like["inobs"]    = "cntcube.fits"
    like["inmodel"]  = "model_bkg.xml"
    like["outmodel"] = "fit_results.xml"
    like["expcube"]  = "expcube.fits"
    like["psfcube"]  = "psfcube.fits"
    like["bkgcube"]  = "bkgcube.fits"
    like["caldb"]    = caldb
    like["irf"]      = irf
    like["debug"]    = True # Switch this always on for results in console
    like.execute()

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Dump header
    print("********************************************")
    print("*      CTA stacked analysis pipeline       *")
    print("********************************************")

    # Setup observations
    obs = setup_observations()

    # Setup model
    obs = setup_model(obs)

    # Run analysis pipeline
    run_pipeline(obs)
