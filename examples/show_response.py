#! /usr/bin/env python
# ==========================================================================
# This script displays one or several CTA response functions.
#
# Required 3rd party modules:
# - matplotlib
# - numpy
#
# Copyright (C) 2014-2015 Juergen Knoedlseder
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
import math
import gammalib
import matplotlib.pyplot as plt


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


# ======================= #
# Show one effective area #
# ======================= #
def show_one_effective_area(rsp, name, color="r"):
    """
    Show one effective area.
    """
    # Set to figure 1
    plt.figure(1)

    # Generate logE(TeV) vector
    nbins = 20
    logE  = [-1.7+0.2*i for i in range(nbins)]
    E     = [math.pow(10.0, e) for e in logE]

    # Generate effective area vector
    aeff = [0.0 for i in range(nbins)]
    for i in range(nbins):
        aeff[i] = rsp.aeff()(logE[i])

    # Plot data
    plt.loglog(E, aeff, color+'-', label=name)
    plt.loglog(E, aeff, color+'o')

    # Set axes
    plt.xlabel("Energy [TeV]")
    plt.ylabel("Effective area [cm2]")

    # Set legend
    plt.legend(loc="lower right")

    # Return
    return


# ======================== #
# Show one background rate #
# ======================== #
def show_one_background_rate(rsp, name, color="r"):
    """
    Show one background rate.
    """
    # Set to figure 2
    plt.figure(2)

    # Generate logE(TeV) vector
    nbins = 20
    logE  = [-1.7+0.2*i for i in range(nbins)]
    E     = [math.pow(10.0, e) for e in logE]

    # Generate on-axis background rate vector
    bgrate = [0.0 for i in range(nbins)]
    for i in range(nbins):
        bgrate[i] = rsp.background()(logE[i], 0.0, 0.0)
        if bgrate[i] < 1.0e-10:
            bgrate[i] = 0.0

    # Plot data
    plt.loglog(E, bgrate, color+'-', label=name)
    plt.loglog(E, bgrate, color+'o')

    # Set axes
    plt.xlabel("Energy [TeV]")
    plt.ylabel("Background rate [events/s/MeV/sr]")

    # Set legend
    plt.legend(loc="lower left")

    # Return
    return


# ==================== #
# Show one sensitivity #
# ==================== #
def show_one_sensitivity(rsp, name, color="r", duration=180000.0, alpha=0.2, sigma=5.0):
    """
    Show one sensitivity.
    """
    # Set constants
    r68_to_sigma = 0.6624305
    TeV2erg      = 1.0e6 * gammalib.MeV2erg

    # Set to figure 3
    plt.figure(3)

    # Generate logE(TeV) vector
    nbins = 20
    logE  = [-1.7+0.2*i for i in range(nbins)]
    E     = [math.pow(10.0, e) for e in logE]

    # Generate sensitivity vector
    flux = [0.0 for i in range(nbins)]
    for i in range(nbins):
        ewidth     = (math.pow(10.0, logE[i]+0.1) - \
                      math.pow(10.0, logE[i]-0.1))*1.0e6
        aeff       = rsp.aeff()(logE[i])
        r68        = rsp.psf().delta_max(logE[i]) / 5.0 / r68_to_sigma
        solidangle = gammalib.twopi * (1.0 - math.cos(r68))
        bgd_counts = rsp.background()(logE[i], 0.0, 0.0) * duration * ewidth * solidangle
        Noff       = bgd_counts/alpha
        if Noff > 0:
            Non        = Non_lima(sigma, Noff)
            src_counts = Non - bgd_counts
            if src_counts < 0.05*bgd_counts:
                src_counts = 0.05*bgd_counts
            if src_counts < 10:
                src_counts = 10.0
            emin    = gammalib.GEnergy(math.pow(10.0, logE[i]-0.1), "TeV")
            emax    = gammalib.GEnergy(math.pow(10.0, logE[i]+0.1), "TeV")
            epivot  = gammalib.GEnergy(math.pow(10.0, logE[i]), "TeV")
            plaw    = gammalib.GModelSpectralPlaw(1.0e-6, -2.6, epivot)
            conv    = TeV2erg*E[i]*E[i]/plaw.flux(emin, emax)
            flux[i] = conv * src_counts / (duration * aeff*0.68)

    # Plot data
    plt.loglog(E, flux, color+'-', label=name)
    plt.loglog(E, flux, color+'o')

    # Set axes
    plt.xlabel("Energy [TeV]")
    plt.ylabel("Sensitivity [erg/cm2/s]")

    # Set legend
    plt.legend(loc="upper right")

    # Return
    return


# ==================== #
# Show one sensitivity #
# ==================== #
def show_one_response(rspname, dbname, name, rootdir=None, color="r"):
    """
    Show one response.
    """
    # Set-up calibration database
    caldb = gammalib.GCaldb()
    if rootdir != None:
        caldb.rootdir(rootdir)
    if gammalib.dir_exists(dbname):
        caldb.rootdir(dbname)
    else:
        caldb.open("cta", dbname)

    # Load response function
    rsp = gammalib.GCTAResponseIrf(rspname, caldb)

    # Show effective area
    show_one_effective_area(rsp, name, color=color)

    # Show background rate
    show_one_background_rate(rsp, name, color=color)

    # Show sensitivity
    show_one_sensitivity(rsp, name, color=color)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create figures
    plt.figure(1)
    plt.title("Effective area")
    plt.figure(2)
    plt.title("Background rate")
    plt.figure(3)
    plt.title("Sensitivity")

    # Set response dictionary
    rsps = [
            {'dbname':  "dummy",
             'rspname': "cta_dummy_irf",
             'name':    "Prod1 (Array E, MPIK)",
             'color':   "b"},
            {'dbname':  "e",
             'rspname': "IFAE20120510_50h",
             'name':    "Prod1 (Array E, IFAE)",
             'color':   "g"},
            {'dbname':  "aar",
             'rspname': "DESY20140105_50h",
             'name':    "Prod2 (Aar, DESY)",
             'color':   "r"}
            ]

    # Loop over all responses
    for rsp in rsps:

        # Get attributes
        dbname  = rsp['dbname']
        rspname = rsp['rspname']
        color   = rsp['color']
        name    = rsp['name']
        if rsp.has_key('rootdir'):
            rootdir = rsp['rootdir']
        else:
            rootdir = None

        # Show response
        show_one_response(rspname, dbname, name, rootdir=rootdir, color=color)

    # Show plot
    plt.show()
