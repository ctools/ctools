#! /usr/bin/env python
# ==========================================================================
# This script peforms a benchmark of the various CTA analysis types
# (unbinned, binned, stacked), and if matplotlib is installed, creates
# a plot of the benchmark results.
#
# Copyright (C) 2014-2016 Jurgen Knodlseder
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
import time
try:
    import matplotlib.pyplot as plt
    has_matplotlib = True
except:
    has_matplotlib = False


# ========================== #
# Unbinned analysis pipeline #
# ========================== #
def unbinned_pipeline(model_name, duration):
    """
    Unbinned analysis pipeline.
    """
    # Set script parameters
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
    cpu_start = time.clock()

    # Simulate events
    sim = ctools.ctobssim()
    sim["inmodel"] = model_name
    sim["caldb"]   = caldb
    sim["irf"]     = irf
    sim["ra"]      = ra
    sim["dec"]     = dec
    sim["rad"]     = rad_sim
    sim["tmin"]    = tstart
    sim["tmax"]    = tstop
    sim["emin"]    = emin
    sim["emax"]    = emax
    sim.run()

    # Select events
    select = ctools.ctselect(sim.obs())
    select["ra"]   = ra
    select["dec"]  = dec
    select["rad"]  = rad_select
    select["tmin"] = tstart
    select["tmax"] = tstop
    select["emin"] = emin
    select["emax"] = emax
    select.run()

    # Get ctlike start CPU time
    cpu_ctlike = time.clock()

    # Perform maximum likelihood fitting
    like = ctools.ctlike(select.obs())
    like.run()

    # Get stop CPU time and compute elapsed times
    cpu_stop    = time.clock()
    cpu_elapsed = cpu_stop - cpu_start
    cpu_ctlike  = cpu_stop - cpu_ctlike

    # Return
    return cpu_elapsed, cpu_ctlike


# ======================== #
# Binned analysis pipeline #
# ======================== #
def binned_pipeline(model_name, duration):
    """
    Binned analysis pipeline.
    """
    # Set script parameters
    caldb       = "prod2"
    irf         = "South_50h"
    ra          =   83.63
    dec         =   22.01
    rad_sim     =   10.0
    tstart      =    0.0
    tstop       = duration
    emin        =    0.1
    emax        =  100.0
    enumbins    =   40
    nxpix       =  200
    nypix       =  200
    binsz       =    0.02
    coordsys    = "CEL"
    proj        = "CAR"

    # Get start CPU time
    cpu_start = time.clock()

    # Simulate events
    sim = ctools.ctobssim()
    sim["inmodel"] = model_name
    sim["caldb"]   = caldb
    sim["irf"]     = irf
    sim["ra"]      = ra
    sim["dec"]     = dec
    sim["rad"]     = rad_sim
    sim["tmin"]    = tstart
    sim["tmax"]    = tstop
    sim["emin"]    = emin
    sim["emax"]    = emax
    sim.run()

    # Bin events into counts map
    bin = ctools.ctbin(sim.obs())
    bin["ebinalg"]  = "LOG"
    bin["emin"]     = emin
    bin["emax"]     = emax
    bin["enumbins"] = enumbins
    bin["nxpix"]    = nxpix
    bin["nypix"]    = nypix
    bin["binsz"]    = binsz
    bin["coordsys"] = coordsys
    bin["xref"]     = ra
    bin["yref"]     = dec
    bin["proj"]     = proj
    bin.run()

    # Get ctlike start CPU time
    cpu_ctlike = time.clock()

    # Perform maximum likelihood fitting
    like = ctools.ctlike(bin.obs())
    like.run()

    # Get stop CPU time and compute elapsed times
    cpu_stop    = time.clock()
    cpu_elapsed = cpu_stop - cpu_start
    cpu_ctlike  = cpu_stop - cpu_ctlike

    # Return
    return cpu_elapsed, cpu_ctlike


# ========================= #
# Stacked analysis pipeline #
# ========================= #
def stacked_pipeline(model_name, duration):
    """
    Stacked analysis pipeline.
    """
    # Set script parameters
    caldb       = "prod2"
    irf         = "South_50h"
    ra          =   83.63
    dec         =   22.01
    rad_sim     =   10.0
    tstart      =    0.0
    tstop       = duration
    emin        =    0.1
    emax        =  100.0
    enumbins    =   40
    nxpix       =  200
    nypix       =  200
    binsz       =    0.02
    coordsys    = "CEL"
    proj        = "CAR"

    # Get start CPU time
    cpu_start = time.clock()

    # Simulate events
    sim = ctools.ctobssim()
    sim["inmodel"] = model_name
    sim["caldb"]   = caldb
    sim["irf"]     = irf
    sim["ra"]      = ra
    sim["dec"]     = dec
    sim["rad"]     = rad_sim
    sim["tmin"]    = tstart
    sim["tmax"]    = tstop
    sim["emin"]    = emin
    sim["emax"]    = emax
    sim.run()

    # Bin events into counts map
    bin = ctools.ctbin(sim.obs())
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
    bin.run()

    # Create exposure cube
    expcube = ctools.ctexpcube(sim.obs())
    expcube["incube"]   = "NONE"
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
    expcube.run()

    # Create PSF cube
    psfcube = ctools.ctpsfcube(sim.obs())
    psfcube["incube"]   = "NONE"
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
    psfcube.run()

    # Create background cube
    bkgcube = ctools.ctbkgcube(sim.obs())
    bkgcube["incube"]   = "NONE"
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
    bkgcube.run()

    # Attach background model to observation container
    bin.obs().models(bkgcube.models())

    # Set Exposure and Psf cube for first CTA observation
    # (ctbin will create an observation with a single container)
    bin.obs()[0].response(expcube.expcube(), psfcube.psfcube(), bkgcube.bkgcube())

    # Get ctlike start CPU time
    cpu_ctlike = time.clock()

    # Perform maximum likelihood fitting
    like = ctools.ctlike(bin.obs())
    like.run()

    # Get stop CPU time and compute elapsed times
    cpu_stop    = time.clock()
    cpu_elapsed = cpu_stop - cpu_start
    cpu_ctlike  = cpu_stop - cpu_ctlike

    # Return
    return cpu_elapsed, cpu_ctlike


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':

    # Initialise parameters
    title      = "CTA analysis benchmark"
    model_name = "${CTOOLS}/share/models/crab.xml"

    # Log header
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
        t_unbinned, ct_unbinned = unbinned_pipeline(model_name, duration)
        t_binned,   ct_binned   = binned_pipeline(model_name, duration)
        t_stacked,  ct_stacked  = stacked_pipeline(model_name, duration)

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
        plt.title(title)
        plt.loglog(a_duration, a_unbinned, 'ro-', label='unbinned')
        plt.loglog(a_duration, a_binned, 'bo-', label='binned')
        plt.loglog(a_duration, a_stacked, 'go-', label='stacked')
        plt.loglog(a_duration, c_unbinned, 'ro--')
        plt.loglog(a_duration, c_binned, 'bo--')
        plt.loglog(a_duration, c_stacked, 'go--')
        plt.xlabel("Observation duration (hours)")
        plt.ylabel("CPU time (seconds)")
        plt.legend(loc="lower right")
        plt.show()
