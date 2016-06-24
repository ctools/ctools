#! /usr/bin/env python
# ==========================================================================
# This script simulates a galactic plane survey of CTA.
#
# Usage:
#   ./make_survey.py
#
# ==========================================================================
#
# Copyright (C) 2011-2015 Juergen Knoedlseder
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
import os
import glob
import sys
import math
import gammalib
from cscripts import obsutils


# =========== #
# Plot counts #
# =========== #
def plot_counts(obs):
    """
    Plot counts.
    """
    # Only proceed if matplotlib is available
    try:
        # Import matplotlib
        import matplotlib.pyplot as plt

        # Set legend fontsize
        params = {'legend.fontsize': 10}
        plt.rcParams.update(params)

        # Set plot styles
        styles = ['b-', 'g-', 'y-', 'n-']

        # Dump header
        print("")
        print("Make plots (using matplotlib):")
        print("==============================")

        # Create figure 1
        plt.figure(1,figsize=(12,6))
        plt.subplots_adjust(hspace=.7)

        # Create subplot 1
        plt.subplot(121)
        plt.title("Spectrum (summed over all pixels)")

        # Loop over observations
        for run in obs:

            # Get event list
            list = run.events()

            # Create energy axis
            ebounds = gammalib.GEbounds()
            emin    = gammalib.GEnergy(0.1, "TeV")
            emax    = gammalib.GEnergy(100.0, "TeV")
            #emin.TeV(0.1)
            #emax.TeV(100.0)
            ebounds.setlog(emin, emax, 10)
            energy = [ebounds.elogmean(i).TeV() for i in range(ebounds.size())]

            # Create spectrum
            print("Extract data:")
            counts = [0.0 for i in range(ebounds.size())]
            print(list.size())
            for atom in list:
                index         = ebounds.index(atom.energy())
                counts[index] = counts[index] + 1.0

            # Create error bars
            error = [math.sqrt(c) for c in counts]

            # Plot spectrum
            plt.loglog(energy, counts, 'ro', label='data')
            #plt.errorbar(energy, counts, error, fmt=None, ecolor='r')

            # Extract models
            print("Extract models:")
            sum_model = [0.0 for i in range(ebounds.size())]
            for k, m in enumerate(obs.models()):
                print("- "+m.name())
                model = [0.0 for i in range(ebounds.size())]
                for atom in list:
                    index        = ebounds.index(atom.energy())
                    prob         = m.eval(atom, run)
                    model[index] = model[index] + prob * atom.size()
                for i in range(ebounds.size()):
                    sum_model[i] = sum_model[i] + model[i]
                #plt.loglog(energy, model, styles[k], label=m.name())
            #plt.loglog(energy, sum_model, 'r-', label='total')

        # Put labels
        plt.xlabel("Energy (TeV)")
        plt.ylabel("Counts")
        plt.legend(loc="lower left")

        # Create subplot 2
        plt.subplot(122)
        plt.title("Offset (summed over all energies)")

        # Set Crab direction
        crab = gammalib.GSkyDir()
        crab.radec_deg(83.63, 22.01)

        # Loop over observations
        for run in obs:

            # Get event list
            list = run.events()

            # Create offset histogram
            print("Extract data:")
            nx       = 30
            doffset2 = 0.01
            offset2  = [(i+0.5)*doffset2 for i in range(nx)]
            counts   = [0.0  for i in range(nx)]
            for atom in list:
                off   = atom.dir().dist_deg(crab)
                off2  = off*off
                index = int(off2/doffset2)
                if index < nx:
                    counts[index] = counts[index] + 1.0

            # Create error bars
            error = [math.sqrt(c) for c in counts]

            # Plot distribution
            #plt.semilogy(offset2, counts, 'ro', label='data')

            # Extract models
            print("Extract models:")
            sum_model = [0.0 for i in range(nx)]
            for k, m in enumerate(obs.models()):
                print("- "+m.name())
                model = [0.0 for i in range(nx)]
                for atom in list:
                    off   = atom.dir().dist_deg(crab)
                    off2  = off*off
                    index = int(off2/doffset2)
                    if index < nx:
                        prob         = m.eval(atom, run)
                        model[index] = model[index] + prob * atom.size()
                for i in range(nx):
                    sum_model[i] = sum_model[i] + model[i]
                #plt.plot(offset2, model, styles[k], label=m.name())
            #plt.plot(offset2, sum_model, 'r-', label='total')
            #plt.ylim(ymin=0.1)

        # Put labels
        plt.xlabel("Offset (deg^2)")
        plt.ylabel("Counts")
        #plt.legend(loc="upper right")

        # Show counts spectra
        plt.show()

    except ImportError:
        print("Matplotlib is not (correctly) installed on your system.")

    # Return
    return


# ======================================== #
# Add CTA background model to observations #
# ======================================== #
def add_background_model(obs):
    """
    Add standard CTA background model to observations container.

    We use a simple power law here, scaled to Konrad's E configuration
    performance table. The model needs still to be validated.
    """
    # Recover models from observation
    models = obs.models()

    # Define background model
    bgd_radial   = gammalib.GCTAModelRadialGauss(3.0)
    bgd_spectrum = gammalib.GModelSpectralPlaw(61.8e-6, -1.85, gammalib.GEnergy(1.0, "TeV"))
    bgd_model    = gammalib.GCTAModelRadialAcceptance(bgd_radial, bgd_spectrum)
    bgd_model.name("Background")
    bgd_model.instruments("CTA")

    # Add background model to container
    models.append(bgd_model)

    # Put container back in observation container
    obs.models(models)

    # Return observation container
    return obs


# ================= #
# Set Crab spectrum #
# ================= #
def crab_spec():
    """
    Set Crab spectrum based on MAGIC observations
    (Albert et al. 2008, ApJ, 674, 1037)
    """
    # Set parameters
    spectrum = gammalib.GModelSpectralPlaw(5.7e-16, -2.48, gammalib.GEnergy(0.3, "TeV"))

    # Return spectrum
    return spectrum


# =============================== #
# Setup single observation survey #
# =============================== #
def survey_single():
    """
    Creates a single observation survey for test purposes.
    """
    # Allocate observation container
    obs = gammalib.GObservations()

    # Set single pointing at galactic centre
    pntdir = gammalib.GSkyDir()
    pntdir.lb_deg(0.0, 0.0)
    run = obsutils.set_obs(pntdir)
    obs.append(run)

    # Define single point source with Crab flux at galactic centre
    center = gammalib.GSkyDir()
    center.lb_deg(0.0, 0.0)
    point_spatial  = gammalib.GModelSpatialPointSource(center)
    point_spectrum = crab_spec()
    point          = gammalib.GModelSky(point_spatial, point_spectrum)
    point.name('GC source')

    # Create model container
    models = gammalib.GModels()
    models.append(point)
    obs.models(models)

    # Return observation container
    return obs


# =========================== #
# Setup Galactic plane survey #
# =========================== #
def survey_gplane(lrange=10, lstep=2):
    """
    Creates a single observation survey for test purposes.

    Keywords:
     lrange - Longitude range (integer deg)
     lstep  - Longitude step size (integer deg)
    """
    # Allocate observation container
    obs = gammalib.GObservations()

    # Loop over longitudes
    for l in range(-lrange,lrange+lstep,lstep):

        # Set pointing
        pntdir = gammalib.GSkyDir()
        pntdir.lb_deg(l, 0.0)
        run = obsutils.set_obs(pntdir)
        run.id(str(l))
        obs.append(run)

    # Define single point source with Crab flux at galactic centre
    center = gammalib.GSkyDir()
    center.lb_deg(0.0, 0.0)
    point_spatial  = gammalib.GModelSpatialPointSource(center)
    point_spectrum = crab_spec()
    point          = gammalib.GModelSky(point_spatial, point_spectrum)
    point.name('GC source')

    # Create model container
    models = gammalib.GModels()
    models.append(point)
    obs.models(models)

    # Return observation container
    return obs


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':

    # Initialise flags
    need_help = False

    # Test for command line arguments
    print(sys.argv[0])
    if (len(sys.argv) > 1):
        if sys.argv[1] == "-h":
            need_help = True
        else:
            need_help = True

    # Print help if needed and exit
    if need_help:
        print("Usage: example_survey.py [OPTIONS]")
        print("     -h       Display this usage message")
        sys.exit()

    # Dump header
    print("***********************")
    print("* Simulate CTA survey *")
    print("***********************")

    # Remove any existing result files
    list = [glob.glob("*.fits"), glob.glob("*.log"), glob.glob("*.xml")]
    for files in list:
        for file in files:
            os.remove(file)

    # Setup single observation survey
    #obs = survey_single()
    obs = survey_gplane()

    # Add background model
    obs = add_background_model(obs)

    # Simulate events
    print("Simulate events")
    obs = obsutils.sim(obs)

    # Make counts map
    print("Make counts map")
    cntmap = obsutils.cntmap(obs)

    # Fit observations
    print("Fit observations")
    like = obsutils.fit(obs)
    #print like.opt()
    #print like.obs().models()

    # Make model map
    print("Make model map (this step will take some time)")
    modmap = obsutils.modmap(obs)

    # Show fit results
    #plot_counts(like.obs())
