#! /usr/bin/env python
# ==========================================================================
# This script peforms a benchmark of the various CTA analysis types
# (unbinned, binned, stacked), and if matplotlib is installed, creates
# a plot of the benchmark results.
#
# Copyright (C) 2014-2015 Jurgen Knodlseder
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
import time
try:
    import matplotlib.pyplot as plt
    has_matplotlib = True
except:
    has_matplotlib = False


# ========================== #
# Unbinned analysis pipeline #
# ========================== #
def unbinned_pipeline(duration):
    """
    Unbinned analysis pipeline.
    """
    # Set script parameters
    model_name  = "${CTOOLS}/share/models/crab.xml"
    caldb       = "prod2"
    irf         = "South_50h"
    ra          =   83.63
    dec         =   22.01
    rad_sim     =   10.0
    tstart      =    0.0
    tstop       = duration
    emin        =    0.1
    emax        =  100.0
    rad_select  =    3.0

    # Get start CPU time
    tstart = time.clock()

    # Simulate events
    sim = ctools.ctobssim()
    sim["inmodel"].filename(model_name)
    sim["caldb"].string(caldb)
    sim["irf"].string(irf)
    sim["ra"].real(ra)
    sim["dec"].real(dec)
    sim["rad"].real(rad_sim)
    sim["tmin"].real(tstart)
    sim["tmax"].real(tstop)
    sim["emin"].real(emin)
    sim["emax"].real(emax)
    sim.run()

    # Select events
    select = ctools.ctselect(sim.obs())
    select["ra"].real(ra)
    select["dec"].real(dec)
    select["rad"].real(rad_select)
    select["tmin"].real(tstart)
    select["tmax"].real(tstop)
    select["emin"].real(emin)
    select["emax"].real(emax)
    select.run()

    # Get ctlike start CPU time
    tctlike = time.clock()

    # Perform maximum likelihood fitting
    like = ctools.ctlike(select.obs())
    like.run()

    # Get stop CPU time
    tstop    = time.clock()
    telapsed = tstop - tstart
    tctlike  = tstop - tctlike
	
    # Return
    return telapsed, tctlike


# ======================== #
# Binned analysis pipeline #
# ======================== #
def binned_pipeline(duration):
    """
    Binned analysis pipeline.
    """
    # Set script parameters
    model_name  = "${CTOOLS}/share/models/crab.xml"
    caldb       = "prod2"
    irf         = "South_50h"
    ra          =   83.63
    dec         =   22.01
    rad_sim     =   10.0
    tstart      =    0.0
    tstop       = duration
    emin        =    0.1
    emax        =  100.0
    enumbins    =   20
    nxpix       =  200
    nypix       =  200
    binsz       =    0.02
    coordsys    = "CEL"
    proj        = "CAR"

    # Get start CPU time
    tstart = time.clock()

    # Simulate events
    sim = ctools.ctobssim()
    sim["inmodel"].filename(model_name)
    sim["caldb"].string(caldb)
    sim["irf"].string(irf)
    sim["ra"].real(ra)
    sim["dec"].real(dec)
    sim["rad"].real(rad_sim)
    sim["tmin"].real(tstart)
    sim["tmax"].real(tstop)
    sim["emin"].real(emin)
    sim["emax"].real(emax)
    sim.run()

    # Bin events into counts map
    bin = ctools.ctbin(sim.obs())
    bin["ebinalg"].string("LOG")
    bin["emin"].real(emin)
    bin["emax"].real(emax)
    bin["enumbins"].integer(enumbins)
    bin["nxpix"].integer(nxpix)
    bin["nypix"].integer(nypix)
    bin["binsz"].real(binsz)
    bin["coordsys"].string(coordsys)
    bin["xref"].real(ra)
    bin["yref"].real(dec)
    bin["proj"].string(proj)
    bin.run()

    # Get ctlike start CPU time
    tctlike = time.clock()

    # Perform maximum likelihood fitting
    like = ctools.ctlike(bin.obs())
    like.run()

    # Get stop CPU time
    tstop    = time.clock()
    telapsed = tstop - tstart
    tctlike  = tstop - tctlike
	
    # Return
    return telapsed, tctlike


# ========================= #
# Stacked analysis pipeline #
# ========================= #
def stacked_pipeline(duration):
    """
    Cube-style analysis pipeline.
    """
    # Set script parameters
    model_name  = "${CTOOLS}/share/models/crab.xml"
    caldb       = "prod2"
    irf         = "South_50h"
    ra          =   83.63
    dec         =   22.01
    rad_sim     =   10.0
    tstart      =    0.0
    tstop       = duration
    emin        =    0.1
    emax        =  100.0
    enumbins    =   20
    nxpix       =  200
    nypix       =  200
    binsz       =    0.02
    coordsys    = "CEL"
    proj        = "CAR"

    # Get start CPU time
    tstart = time.clock()

    # Simulate events
    sim = ctools.ctobssim()
    sim["inmodel"].filename(model_name)
    sim["caldb"].string(caldb)
    sim["irf"].string(irf)
    sim["ra"].real(ra)
    sim["dec"].real(dec)
    sim["rad"].real(rad_sim)
    sim["tmin"].real(tstart)
    sim["tmax"].real(tstop)
    sim["emin"].real(emin)
    sim["emax"].real(emax)
    sim.run()

    # Bin events into counts map
    bin = ctools.ctbin(sim.obs())
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
    bin.run()

    # Create exposure cube
    expcube = ctools.ctexpcube(sim.obs())
    expcube["incube"].filename("NONE")
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
    expcube.run()

    # Create PSF cube
    psfcube = ctools.ctpsfcube(sim.obs())
    psfcube["incube"].filename("NONE")
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
    psfcube.run()

    # Create background cube
    bkgcube = ctools.ctbkgcube(sim.obs())
    bkgcube["incube"].filename("NONE")
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
    bkgcube.run()

    # Attach background model to observation container
    bin.obs().models(bkgcube.models())

    # Set Exposure and Psf cube for first CTA observation
    # (ctbin will create an observation with a single container)
    bin.obs()[0].response(expcube.expcube(), psfcube.psfcube(), bkgcube.bkgcube())

    # Get ctlike start CPU time
    tctlike = time.clock()

    # Perform maximum likelihood fitting
    like = ctools.ctlike(bin.obs())
    like.run()

    # Get stop CPU time
    tstop    = time.clock()
    telapsed = tstop - tstart
    tctlike  = tstop - tctlike
	
    # Return
    return telapsed, tctlike


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    """
    # Dump header
    print("*************************************")
    print("*      CTA analysis benchmark       *")
    print("*************************************")
    print("  Duration   Unbinned     Binned    Stacked  ctlike-UB   ctlike-B   ctlike-S")

    # Initialize arrays
    a_duration = []
    a_unbinned = []
    a_binned   = []
    a_stacked  = []
    c_unbinned = []
    c_binned   = []
    c_stacked  = []

    # Perform benchmarks
    duration = 1800.0
    for i in range(10):

        # Perform analyses
        t_unbinned, ct_unbinned = unbinned_pipeline(duration)
        t_binned,   ct_binned   = binned_pipeline(duration)
        t_stacked,  ct_stacked  = stacked_pipeline(duration)

        # Collect results
        a_duration.append(duration/3600.0)
        a_unbinned.append(t_unbinned)
        a_binned.append(t_binned)
        a_stacked.append(t_stacked)
        c_unbinned.append(ct_unbinned)
        c_binned.append(ct_binned)
        c_stacked.append(ct_stacked)

        # Print results
        print("%10.0f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f" % \
              (duration, t_unbinned, t_binned, t_stacked, \
                         ct_unbinned, ct_binned, ct_stacked))

        # Increment duration
        duration *= 2.0
 
    # Optionally plot results
    if has_matplotlib:
        plt.figure(1)
        plt.title("CTA analysis benchmark")
        plt.loglog(a_duration, a_unbinned, 'ro-', label='unbinned')
        plt.loglog(a_duration, a_binned, 'bo-', label='binned')
        plt.loglog(a_duration, a_stacked, 'go-', label='cube-style')
        plt.loglog(a_duration, c_unbinned, 'ro--')
        plt.loglog(a_duration, c_binned, 'bo--')
        plt.loglog(a_duration, c_stacked, 'go--')
        plt.xlabel("Duration (hours)")
        plt.ylabel("CPU time (seconds)")
        plt.legend(loc="lower right")
        plt.show()
        