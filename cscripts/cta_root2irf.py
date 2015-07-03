#! /usr/bin/env python
# ==========================================================================
# This script generates ASCII performance files and file functions for
# background modelling from MC ROOT files. On input it takes the filename
# of the ROOT performance file while on output it produced 2 ASCII files.
# The file "performance.dat" contains the performance data in the usual
# ASCII format, the file "background.txt" contains a file function that can
# be used within GammaLib for spectral background modelling.
#
# Note that the effective area given in the ASCII file is the full effective
# area (i.e. the area before applying any theta cut). This is required by
# GammaLib as the maximum likelihood analysis does not imply any theta cut.
#
# Copyright (C) 2011-2014 Juergen Knoedlseder
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
from ROOT import TFile, TH1F, TH2F
import math
import sys
import glob
import os


# ================================== #
# Compute Eq. (17) of Li & Ma (1983) #
# ================================== #
def sigma_lima(Non, Noff, alpha=0.2):
    """
    Compute Eq. (17) of Li & Ma (1983).

    Parameters:
     Non   - Number of on counts
     Noff  - Number of off counts
    Keywords:
     alpha - Ratio of on-to-off exposure
    """
    # Compute sensitivity
    alpha1 = alpha + 1.0
    sum    = Non + Noff
    arg1   = Non / sum
    arg2   = Noff / sum
    term1  = Non  * math.log((alpha1/alpha)*arg1)
    term2  = Noff * math.log(alpha1*arg2)
    sigma  = math.sqrt(2.0 * (term1 + term2))

    # Return sensitivity
    return sigma


# ======================================== #
# Solve Eq. (17) of Li & Ma (1983) for Non #
# ======================================== #
def Non_lima(sigma, Noff, alpha=0.2):
    """
    Solve Eq. (17) of Li & Ma (1983) for Non

    Parameters:
     sigma - Required significance
     Noff  - Number of off counts
    Keywords:
     alpha - Ratio of on-to-off exposure
    """
    # Set initial guess for Non
    S   = 25.0
    Non = S + alpha*Noff

    # Compute initial sigma estimate
    s     = sigma_lima(Non, Noff, alpha=alpha)
    delta = s - sigma

    # Iterate
    while abs(delta) > 1.0e-6:
        S    *= sigma/s
        Non   = S + alpha*Noff
        s     = sigma_lima(Non, Noff, alpha=alpha)
        delta = s - sigma

    # Return
    return Non


# =============================================== #
# Compute effective area assuming a Crab spectrum #
# =============================================== #
def get_aeff(bgrate, emin, emax, \
             duration=180000.0, alpha=0.2, sigma=5.0, norm=2.79e-7, index=-2.57):
    """
    Computes effective area assuming a Crab spectrum. This does not
    include any hadron cut efficiencies, and thus is not very useful in a first
    place ...

    Parameters:
     bgrate   - Background rate (TeV)
     emin     - Minimum energy (TeV)
     emax     - Maximum energy (TeV)
    Keywords:
     duration - Observing time in seconds (default: 50h)
     alpha    - Non/Noff ratio (default: 0.2)
     sigma    - Significance (default: 5.0)
     norm     - Crab spectrum normalization at 1 TeV (default: 2.79e-7)
     index    - Crab spectral index (default: -2.57)
    """
    # Compute number of source events
    B    = bgrate*duration
    Noff = B/alpha
    Non  = Non_lima(sigma, Noff)
    S    = (Non - B)

    # Compute Crab flux in specified energy band
    inx = index+1.0
    cu  = norm/inx * (math.pow(emax, inx) - math.pow(emin, inx))

    # Compute effective area
    area = S / duration / cu

    # Return effective area
    return area


# ============== #
# Generate files #
# ============== #
def make_files(filename, prfname, bkgname):
    """
    Generate files.

    Parameters:
     filename - Performance file name
     prfname  - Performance file name
     bkgname  - File function filename
    """
    # Open performance file and file function
    file = TFile(filename)
    fprf = open(prfname, "w")
    fbkg = open(bkgname, "w")

    # Write header
    fprf.write("log(E) Area  r68  r80  ERes. BG Rate  Diff Sens\n")
    #sys.stdout.write("log(E) Area  Area80  r68  r80  ERes. BG Rate  Diff Sens  BGRateSqDeg  BGInt  (BGIntControl)\n")

    # Get histograms. Konrad's root files do not have an EffectiveArea80
    # histogram, but his effective area is for 80% containment radius
    # anyways, so we simply read this histogram into aeff80.
    sens        = TH1F()
    bgrate      = TH1F()
    bgratesqdeg = TH1F()
    aeff        = TH1F()
    aeff80      = TH1F()
    angres      = TH1F()
    angres80    = TH1F()
    eres        = TH1F()
    file.GetObject("DiffSens", sens)
    file.GetObject("BGRate", bgrate)
    file.GetObject("BGRatePerSqDeg", bgratesqdeg)
    file.GetObject("EffectiveArea", aeff)
    try:
        file.GetObject("EffectiveArea80", aeff80)
    except:
        file.GetObject("EffectiveArea", aeff80)
    file.GetObject("AngRes", angres)
    file.GetObject("AngRes80", angres80)
    file.GetObject("ERes", eres)

    # Setup vectors
    energies  = bgrate.GetXaxis()
    nbins_eng = energies.GetNbins()
    for i in range(nbins_eng):

        # Get logE and background rate per squared degree
        logE   = energies.GetBinCenter(i+1)
        bgd    = bgrate.GetBinContent(i+1,1)
        bgdsq  = bgratesqdeg.GetBinContent(i+1,1)
        area   = aeff.GetBinContent(i+1,1)
        area80 = aeff80.GetBinContent(i+1,1)
        r68    = angres.GetBinContent(i+1,1)
        r80    = angres80.GetBinContent(i+1,1)
        dE     = eres.GetBinContent(i+1,1)
        diffs  = sens.GetBinContent(i+1,1)

        # Skip any NaNs
        if math.isnan(area80):
            continue

        # Compute energy (in MeV)
        energy = math.pow(10.0, logE)*1.0e6
        emin   = math.pow(10.0, logE-0.1)*1.0e6
        emax   = math.pow(10.0, logE+0.1)*1.0e6
        ewidth = emax-emin

        # Convert into background rate per steradian and MeV
        bkg_rate  = bgdsq / (0.01745329*0.01745329) / ewidth

        # Compute control background rate (this only works for KB files
        # as they are for 80% containment radius). But we don't use this
        # anyways, it's just for control on display ...
        omega = 2.0 * math.pi * (1.0 - math.cos(math.radians(r80)))
        if omega > 0.0:
            bkg_rate2 = bgd / omega / ewidth
        else:
            bkg_rate2 = 0.0

        # Compute full effective area by dividing area80/0.80
        area = area80/0.80

        # Write results in file
        line = "%.1f  %.1f  %.4f  %.4f %.4f %.5e  %.5e" % \
              (logE, area, r68, r80, dE, bgd, diffs)
        fprf.write(line+"\n")
        fbkg.write(str(energy)+" "+str(bkg_rate)+"\n")

        # Show performance file
        #print "%.1f  %.1f  %.1f  %.4f  %.4f %.4f %.7f  %.5e  %.5e %.5e (%.5e)" % \
        #      (logE, area, area80, r68, r80, dE, bgd, diffs, bgdsq, bkg_rate, bkg_rate2)

    # Write trailer
    fprf.write("---------------------------------------------\n")
    fprf.write("Notes\n")
    fprf.write(" 1) log(E) = log10(E/TeV) - bin centre\n")
    fprf.write(" 2) Eff Area - in square metres after background cut (no theta cut)\n")
    fprf.write(" 3) Ang. Res - 68% containment radius of gamma-ray PSF post cuts - in degrees\n")
    fprf.write(" 4) Ang. Res - 80% containment radius of gamma-ray PSF post cuts - in degrees\n")
    fprf.write(" 5) Fractional Energy Resolution (rms)\n")
    fprf.write(" 6) BG Rate  - inside point-source selection region - post call cuts - in Hz\n")
    fprf.write(" 7) Diff Sens - differential sensitivity for this bin expressed as E^2 dN/dE\n")
    fprf.write(" - in erg cm^-2 s^-1 - for a 50 hours exposure - 5 sigma significance including\n")
    fprf.write(" systematics and statistics and at least 10 photons.\n")

    # Close files
    fprf.close()
    fbkg.close()

    # Return
    return


#==========================#
# Main routine entry point #
#==========================#
if __name__ == '__main__':
    """
    """
    # Get ROOT filename from command line argument
    if (len(sys.argv) > 1):
        filename = sys.argv[1]
    else:
        #caldb = "/Users/jurgen/Documents/Travail/Projects/CTA/WP-MC/root/IFAEPerformanceBCE_Nov2010"
        #caldb = "/Users/jurgen/Documents/Travail/Projects/CTA/WP-MC/root/Kb_subarrays_root"
        caldb = "/Users/jurgen/Documents/Travail/Projects/CTA/WP-MC/root/Kb_all_50deg_root"

    # Get list of all performance files
    files = glob.glob(caldb+"/*.root")
    files.sort()

    # Loop over all files
    for file in files:

        # Get filename
        head, tail = os.path.split(file)

        # Build name
        name    = tail.strip(".root")
        if name.find("kb_") != -1:
            name = name.replace("_20deg", "")
        elif name.find("IFAE_") != -1:
            name = name.replace("Subarray", "")
            name = name.replace("IFAE_", "")
            name = "IFAE_"+name
            i    = name.find("hours")
            if i != -1:
                name = name[0:i+1]

        # Build output filenames
        prfname = name + ".dat"
        bkgname = "bkg_"+name+".txt"
        #print prfname, bkgname

        # Create files
        make_files(file, prfname, bkgname)
