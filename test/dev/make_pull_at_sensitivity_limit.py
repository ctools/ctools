#! /usr/bin/env python
# ==========================================================================
# This script generates pull distributions for a source at the 5 sigma
# sensitivity limit. The script makes use of the "processing" Python module,
# if available, that can be used for parallel computing on multiple
# cores/CPUs.
#
# Required 3rd party modules:
# - processing (optional)
#
# Copyright (C) 2011 Jurgen Knodlseder
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
import sys
import time
import csv
import gammalib
from cscripts import cspull

# Try importing processing module
has_processing = False
try:
    import processing
    has_processing = True
except:
    pass


# ======================== #
# Create pull distribution #
# ======================== #
def create_pull(loge, emin, emax, models, ntrials=100, duration=180000.0,
                enumbins=0, log=False):
    """
    Create pull distribution

    Parameters:
     emin - Minimum energy (TeV)
     emax - Maximum energy (TeV)
    Keywords:
     log  - Create log file(s)
    """
    # Generate output filename
    outfile = "pull_"+loge+".dat"

    # Setup cspull tool
    pull = cspull()
    pull.models(models)
    pull["outfile"]  = outfile
    pull["ntrials"]  = ntrials
    pull["caldb"]    = "prod2"
    pull["irf"]      = "South_50h"
    pull["ra"]       = 83.6331
    pull["dec"]      = 22.0145
    pull["emin"]     = emin
    pull["emax"]     = emax
    pull["enumbins"] = enumbins
    pull["duration"] = duration

    # Optionally open the log file
    if log:
        pull.logFileOpen()

    # Run tool
    pull.run()

    # Return
    return


# ============ #
# Create model #
# ============ #
def create_model(flux, index=-2.48, fitidx=False):
    """
    Add standard CTA background model to observations container.
    We use a simple power law here, scaled to Konrad's E configuration
    performance table.

    Parameters:
     flux - Flux in Crab units
    """
    # Define background model
    bgd_radial   = gammalib.GCTAModelRadialGauss(3.0)
    bgd_spectrum = gammalib.GModelSpectralPlaw(61.8, -1.85)
    bgd_spectrum["Prefactor"].scale(1.0e-6)
    bgd_spectrum["PivotEnergy"].value(1.0)
    bgd_spectrum["PivotEnergy"].scale(1.0e6)
    if fitidx:
        bgd_spectrum["Index"].free()
    else:
        bgd_spectrum["Index"].fix()
    bgd_model = gammalib.GCTAModelRadialAcceptance(bgd_radial, bgd_spectrum)
    bgd_model.name("Background")
    bgd_model.instruments("CTA")

    # Define source spectrum
    location = gammalib.GSkyDir()
    location.radec_deg(83.6331, 22.0145)
    src_spatial  = gammalib.GModelSpatialPtsrc(location)
    src_spectrum = gammalib.GModelSpectralPlaw(flux, index)
    src_spectrum["Prefactor"].scale(5.7e-16)
    src_spectrum["PivotEnergy"].value(0.3)
    src_spectrum["PivotEnergy"].scale(1.0e6)
    if fitidx:
        src_spectrum["Index"].free()
    else:
        src_spectrum["Index"].fix()	
    src_model = gammalib.GModelPointSource(src_spatial, src_spectrum)
    src_model.name("Test")

    # Add models to container
    models = gammalib.GModels()
    models.append(bgd_model)
    models.append(src_model)

    # Return model container
    return models


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Get input arguments
    usage = """
    make_pull_at_sensitivity_limit filename [max_threads]
    Run cssens to create the input file.
    """

    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print(usage)
        sys.exit()

    # Get sensitivity file
    filename = sys.argv[1]

    # Set maximum number of threads
    if len(sys.argv) > 2:
        max_threads = int(sys.argv[2])
    else:
        max_threads = 1

    # Load sensitivity file
    file = open(filename, 'r')
    reader = csv.DictReader(file)

    # Loop over rows
    for row in reader:

        # Get parameters
        loge = float(row['loge'])
        if loge < 0:
            loge = "m"+str(abs(loge))
        else:
            loge = "p"+str(abs(loge))
        emin = float(row['emin'])
        emax = float(row['emax'])
        flux = float(row['crab'])

        # Generate model
        models = create_model(flux)

        # Processing support?
        if has_processing:

            # Wait until one thread has finished
            while len(processing.activeChildren()) >= max_threads:
                time.sleep(60)

            # Set arguments
            args   = (loge, emin, emax, models)
            kwargs = {}

            # Generate pull distribution
            p = processing.Process(target=create_pull, args=args, kwargs=kwargs)
            p.start()
            print("Process emin=%.4f emax=%.4f started." % (emin, emax))

            # Wait some time
            time.sleep(1)

        # ... no
        else:
            create_pull(loge, emin, emax, models)

    # Processing support
    if has_processing:

        # Wait until all threads finished
        while len(processing.activeChildren()) > 0:
            time.sleep(60)

    # Close file
    file.close()
