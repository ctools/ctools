#! /usr/bin/env python
# ==========================================================================
# This script illustrates how to perform a stacked cube analysis.
#
# Copyright (C) 2014 Juergen Knoedlseder
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
import ctools
import gammalib


# ================= #
# Analysis pipeline #
# ================= #
def pipeline(caldb = "dummy", \
             irf   = "cta_dummy_irf", \
             model = "$CTOOLS/share/models/crab.xml"):
    """
    Analysis pipeline for a stacked cube analysis. The pipeline writes
    intermediate results to FITS files on disk. Also log files of all
    individual ctools will be written.
    
    Keywords:
     caldb - Calibration database repository
     irf   - Name of instrument response
     model - Source model
    """
    # Set pointing pattern (4 pointings)
    pointings = ctools.obsutils.set_obs_patterns("four", offset=2.0)

    # Set CTA observations (default: 1800 sec per pointing)
    obs = ctools.obsutils.set_obs_list(pointings, caldb=caldb, irf=irf)

    # Load models
    models = gammalib.GModels(model)
    obs.models(models)

    # Run ctobssim (simulates 4 event lists)
    sim = ctools.ctobssim(obs)
    sim.logFileOpen()
    sim.run()

    # Run ctbin (builds stacked event cube)
    bin = ctools.ctbin(sim.obs())
    bin.logFileOpen()
    bin["outfile"].filename("cntmap.fits")
    bin["ebinalg"].string("LOG")
    bin["emin"].real(0.1)
    bin["emax"].real(100.0)
    bin["enumbins"].integer(20)
    bin["usepnt"].boolean(False)
    bin["nxpix"].integer(300)
    bin["nypix"].integer(300)
    bin["binsz"].real(0.02)
    bin["coordsys"].string("CEL")
    bin["proj"].string("CAR")
    bin["xref"].real(83.6331)
    bin["yref"].real(22.0145)
    bin.execute()

    # Run ctexpcube (builds exposure cube)
    expcube = ctools.ctexpcube(sim.obs())
    expcube.logFileOpen()
    expcube["cntmap"].filename("cntmap.fits")
    expcube["outfile"].filename("expcube.fits")
    expcube["caldb"].string(caldb)
    expcube["irf"].string(irf)
    expcube.execute()

    # Run ctpsfcube (builds PSF cube)
    psfcube = ctools.ctpsfcube(sim.obs())
    psfcube.logFileOpen()
    psfcube["cntmap"].filename("NONE")
    psfcube["outfile"].filename("psfcube.fits")
    psfcube["caldb"].string(caldb)
    psfcube["irf"].string(irf)
    psfcube["ebinalg"].string("LOG")
    psfcube["emin"].real(0.1)
    psfcube["emax"].real(100.0)
    psfcube["enumbins"].integer(20)
    psfcube["nxpix"].integer(5)
    psfcube["nypix"].integer(5)
    psfcube["binsz"].real(2.0)
    psfcube["coordsys"].string("CEL")
    psfcube["proj"].string("CAR")
    psfcube["xref"].real(83.6331)
    psfcube["yref"].real(22.0145)
    psfcube["amax"].real(0.3)
    psfcube["anumbins"].integer(60)
    psfcube.execute()

    # Run ctbkgcube (builds background cube)
    bkgcube = ctools.ctbkgcube(sim.obs())
    bkgcube.logFileOpen()
    bkgcube["cntmap"].filename("cntmap.fits")
    bkgcube["bkgmdl"].filename(model)
    bkgcube["outfile"].filename("bkgcube.fits")
    bkgcube.execute()

    # Create observation container for stacked cube analysis
    obs.clear()
    cta = gammalib.GCTAObservation("cntmap.fits", "expcube.fits", "psfcube.fits")
    obs.append(cta)

    # Create model for stacked cube analysis
    models = gammalib.GModels(model)
    models.remove("Background")
    spatial    = gammalib.GModelSpatialDiffuseCube("bkgcube.fits")
    spectral   = gammalib.GModelSpectralConst()
    background = gammalib.GCTAModelCubeBackground(spatial, spectral)
    models.append(background)
    models.save("cube_model.xml")

    # Append model to observation container
    obs.models(models)

    # Run ctlike
    like = ctools.ctlike(obs)
    like.logFileOpen()
    like["debug"].boolean(True)
    like["outmdl"].filename("cube_results.xml")
    like.execute()

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    """
    # Dump header
    print("*************************************")
    print("*     CTA stacked cube analysis     *")
    print("*************************************")

    # Run pipeline
    pipeline()
    