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
from ctools import obsutils


# ================== #
# Setup observations #
# ================== #
def setup_observations(pattern="four", ra=83.63, dec=22.01, offset=1.5, \
                       emin=0.1, emax=100.0, rad=5.0, duration=1800.0, \
                       deadc=0.95, \
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
    obs_def_list = obsutils.set_obs_patterns(pattern, \
                                             ra=ra, \
                                             dec=dec, \
                                             offset=offset)
    
    # Get observation container
    obs = obsutils.set_obs_list(obs_def_list, \
                                duration=duration, \
                                emin=emin, \
                                emax=emax, \
                                rad=rad, \
                                caldb=caldb, \
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
def run_pipeline(obs, ra=83.63, dec=22.01, emin=0.1, emax=100.0, \
                 enumbins=20, nxpix=200, nypix=200, binsz=0.02, \
                 coordsys="CEL", proj="CAR", \
                 model="${CTOOLS}/share/models/crab.xml", \
                 caldb="prod2", irf="South_50h", \
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
    sim["debug"].boolean(debug)
    sim["outevents"].filename("obs.xml")
    sim.execute()

    # Bin events into counts map
    bin = ctools.ctbin()
    bin["inobs"].filename("obs.xml")
    bin["outcube"].filename("cntcube.fits")
    bin["ebinalg"].string("LOG")
    bin["emin"].real(emin)
    bin["emax"].real(emax)
    bin["enumbins"].integer(enumbins)
    bin["nxpix"].integer(nxpix)
    bin["nypix"].integer(nypix)
    bin["binsz"].real(binsz)
    bin["coordsys"].string(coordsys)
    bin["proj"].string(proj)
    bin["xref"].real(ra)
    bin["yref"].real(dec)
    bin["debug"].boolean(debug)
    bin.execute()

    # Create exposure cube
    expcube = ctools.ctexpcube()
    expcube["inobs"].filename("obs.xml")
    expcube["incube"].filename("cntcube.fits")
    expcube["outcube"].filename("expcube.fits")
    expcube["caldb"].string(caldb)
    expcube["irf"].string(irf)
    expcube["ebinalg"].string("LOG")
    expcube["emin"].real(emin)
    expcube["emax"].real(emax)
    expcube["enumbins"].integer(enumbins)
    expcube["nxpix"].integer(nxpix)
    expcube["nypix"].integer(nypix)
    expcube["binsz"].real(binsz)
    expcube["coordsys"].string(coordsys)
    expcube["proj"].string(proj)
    expcube["xref"].real(ra)
    expcube["yref"].real(dec)
    expcube["debug"].boolean(debug)
    expcube.execute()

    # Create PSF cube
    psfcube = ctools.ctpsfcube()
    psfcube["inobs"].filename("obs.xml")
    psfcube["incube"].filename("NONE")
    psfcube["outcube"].filename("psfcube.fits")
    psfcube["caldb"].string(caldb)
    psfcube["irf"].string(irf)
    psfcube["ebinalg"].string("LOG")
    psfcube["emin"].real(emin)
    psfcube["emax"].real(emax)
    psfcube["enumbins"].integer(enumbins)
    psfcube["nxpix"].integer(10)
    psfcube["nypix"].integer(10)
    psfcube["binsz"].real(1.0)
    psfcube["coordsys"].string(coordsys)
    psfcube["proj"].string(proj)
    psfcube["xref"].real(ra)
    psfcube["yref"].real(dec)
    psfcube["debug"].boolean(debug)
    psfcube.execute()

    # Create background cube
    bkgcube = ctools.ctbkgcube()
    bkgcube["inobs"].filename("obs.xml")
    bkgcube["inmodel"].filename(model)
    bkgcube["incube"].filename("cntcube.fits")
    bkgcube["outcube"].filename("bkgcube.fits")
    bkgcube["outmodel"].filename("model_bkg.xml")
    bkgcube["caldb"].string(caldb)
    bkgcube["irf"].string(irf)
    bkgcube["ebinalg"].string("LOG")
    bkgcube["emin"].real(emin)
    bkgcube["emax"].real(emax)
    bkgcube["enumbins"].integer(enumbins)
    bkgcube["nxpix"].integer(10)
    bkgcube["nypix"].integer(10)
    bkgcube["binsz"].real(1.0)
    bkgcube["coordsys"].string(coordsys)
    bkgcube["proj"].string(proj)
    bkgcube["xref"].real(ra)
    bkgcube["yref"].real(dec)
    bkgcube["debug"].boolean(debug)
    bkgcube.execute()

    # Perform maximum likelihood fitting
    like = ctools.ctlike()
    like["inobs"].filename("cntcube.fits")
    like["inmodel"].filename("model_bkg.xml")
    like["outmodel"].filename("fit_results.xml")
    like["expcube"].filename("expcube.fits")
    like["psfcube"].filename("psfcube.fits")
    like["bkgcube"].filename("bkgcube.fits")
    like["caldb"].string(caldb)
    like["irf"].string(irf)
    like["debug"].boolean(True) # Switch this always on for results in console
    like.execute()
	
    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    """
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
