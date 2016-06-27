#! /usr/bin/env python
# ==========================================================================
# This script simulates events from a data-space model and compares them to
# the model.
#
# If matplotlib is installed, the event spectrum will be displayed.
#
# --------------------------------------------------------------------------
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
import gammalib
import ctools
import math
try:
    import matplotlib.pyplot as plt
    has_matplotlib = True
except:
    has_matplotlib = False


# ======================== #
# Create CTA observation #
# ======================== #
def createobs(ra=86.171648, dec=-1.4774586, rad=5.0,
              emin=0.1, emax=100.0, duration=360000.0, deadc=0.95,
              irf="South_50h", caldb="prod2"):
    """
    Create CTA observation.
    """
    # Allocate CTA observation
    obs = gammalib.GCTAObservation()

    # Set calibration database
    db = gammalib.GCaldb()
    if (gammalib.dir_exists(caldb)):
        db.rootdir(caldb)
    else:
        db.open("cta", caldb)

    # Set pointing direction
    pntdir = gammalib.GSkyDir()
    pntdir.radec_deg(ra, dec)
    pnt = gammalib.GCTAPointing()
    pnt.dir(pntdir)
    obs.pointing(pnt)

    # Set ROI
    roi     = gammalib.GCTARoi()
    instdir = gammalib.GCTAInstDir()
    instdir.dir(pntdir)
    roi.centre(instdir)
    roi.radius(rad)

    # Set GTI
    gti   = gammalib.GGti()
    start = gammalib.GTime(0.0)
    stop  = gammalib.GTime(duration)
    gti.append(start, stop)

    # Set energy boundaries
    ebounds = gammalib.GEbounds()
    e_min   = gammalib.GEnergy()
    e_max   = gammalib.GEnergy()
    e_min.TeV(emin)
    e_max.TeV(emax)
    ebounds.append(e_min, e_max)

    # Allocate event list
    events = gammalib.GCTAEventList()
    events.roi(roi)
    events.gti(gti)
    events.ebounds(ebounds)
    obs.events(events)

    # Set instrument response
    obs.response(irf, db)

    # Set ontime, livetime, and deadtime correction factor
    obs.ontime(duration)
    obs.livetime(duration*deadc)
    obs.deadc(deadc)

    # Return observation
    return obs    


# =============== #
# Simulate events #
# =============== #
def simulate(obs, xmlname):
    """
    Simulate events.
    """
    # Allocate random number generator
    ran = gammalib.GRan()

    # Load models and extract first model
    models = gammalib.GModels(xmlname)
    model = models[0]
    print(model)

    # Simulate events
    events = model.mc(obs, ran)

    # Print event statistics
    npred = obs.npred(models)
    print(str(len(events)) + " events simulated.")
    print(str(npred) + " events expected from Npred.")

    # Check distance
    pntdir = obs.pointing().dir()
    for event in events:
        distance = gammalib.GCTAInstDir(event.dir()).dir().dist_deg(pntdir)
        if (distance > 5.0):
            print("Too far: "+str(distance))

    # Return events
    return events


# ============================== #
# Simulate events using ctobssim #
# ============================== #
def simulate_ctobssim(obs, xmlname, seed=0):
    """
    Simulate events using ctobssim.
    """
    # Create containers
    observations = gammalib.GObservations()
    observations.append(obs)

    # Append models
    models = gammalib.GModels(xmlname)
    observations.models(models)
    print(models[0])

    # Allocate ctobssim application and set parameters
    sim = ctools.ctobssim(observations)
    sim['seed'].integer(seed)

    # Run simulator
    sim.run()

    # Retrieve events
    events = sim.obs()[0].events().copy()

    # Print event statistics
    npred = obs.npred(models)
    print(str(len(events)) + " events simulated.")
    print(str(npred) + " events expected from Npred.")

    # Delete the simulation
    del sim

    # Return events
    return events


# =========== #
# Show events #
# =========== #
def show_events(events, xmlname, duration, emin, emax, ebins=30):
    """
    Show events using matplotlib.
    """
    # Create figure
    plt.figure(1)
    plt.title("MC simulated event spectrum (" + str(emin) + '-' + str(emax) + " TeV)")

    # Setup energy range covered by data
    ebds = gammalib.GEbounds(ebins, gammalib.GEnergy(emin, "TeV"),
                                    gammalib.GEnergy(emax, "TeV"))

    # Create energy axis
    energy = []
    for i in range(ebds.size()):
        energy.append(ebds.elogmean(i).TeV())

    # Fill histogram
    counts = [0.0 for i in range(ebds.size())]
    for event in events:
        index = ebds.index(event.energy())
        counts[index] = counts[index] + 1.0

    # Create error bars
    error = [math.sqrt(c) for c in counts]

    # Get model values
    sum    = 0.0
    models = gammalib.GModels(xmlname)
    m      = models[0]
    model  = []
    t = gammalib.GTime()
    for i in range(ebds.size()):
        eval   = ebds.elogmean(i)
        ewidth = ebds.emax(i) - ebds.emin(i)
        f      = m.npred(eval, t, obs) * ewidth.MeV() * duration
        sum   += f
        model.append(f)
    print(str(sum) + " events expected from spectrum (integration).")

    # Plot data
    plt.loglog(energy, counts, 'ro')
    plt.errorbar(energy, counts, error, fmt=None, ecolor='r')

    # Plot model
    plt.plot(energy, model, 'b-')

    # Set axes
    plt.xlabel("Energy (TeV)")
    plt.ylabel("Number of events")

    # Notify
    print("PLEASE CLOSE WINDOW TO CONTINUE ...")

    # Show plot
    plt.show()

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':

    # Dump header
    print('')
    print('*******************')
    print('* Simulate events *')
    print('*******************')

    # Set XML names
    xmlname = 'model_bgd.xml'

    # Set simulation parameters
    duration = 360000.0 # Observation duration (s)
    emin     = 0.1      # 0.1 TeV
    emax     = 100.0    # 100 TeV
    ebins    = 30

    # Set observation
    obs = createobs(emin=emin, emax=emax, duration=duration)

    # Perform simulation
    #events = simulate(obs, xmlname)
    events = simulate_ctobssim(obs, xmlname)

    # Show events
    if has_matplotlib:
        show_events(events, xmlname, duration, emin, emax, ebins=ebins)
