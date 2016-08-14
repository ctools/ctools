#! /usr/bin/env python
# ==========================================================================
# Shows one or several CTA response functions.
#
# Copyright (C) 2014-2016 Juergen Knoedlseder
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
import cscripts
try:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.close()
except (ImportError, RuntimeError):
    print('This script needs the "matplotlib" module')
    sys.exit()


# ================================== #
# Compute Eq. (17) of Li & Ma (1983) #
# ================================== #
def sigma_lima(Non, Noff, alpha=0.2):
    """
    Compute Eq. (17) of Li & Ma (1983).

    Parameters
    ----------
    Non : float
        Number of On counts
    Noff : float
        Number of Off counts
    alpha : float, optional
        Ratio of on-to-off exposure

    Returns
    -------
    sigma : float
        Detection significance in Gaussian sigma
    """
    # Compute sensitivity in Gaussian sigma
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

    Parameters
    ----------
    sigma : float
        Required detection significance in Gaussian sigma
    Noff : float
        Number of Off counts
    alpha : float, optional
        Ratio of on-to-off exposure

    Returns
    -------
    Non : float
        Number of On counts
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
def show_one_effective_area(rsp, name, color='r'):
    """
    Show one effective area

    Parameters
    ----------
    rsp : `~gammalib.GCTAResponse`
        Response function
    name : str
        Name of the response function
    color : str, optional
        Color for plot
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
    plt.xlabel('Energy (TeV)')
    plt.ylabel(r'Effective area (cm$^{2}$)')

    # Set legend
    plt.legend(loc='lower right')

    # Return
    return


# ======================== #
# Show one background rate #
# ======================== #
def show_one_background_rate(rsp, name, color='r'):
    """
    Show one background rate

    Parameters
    ----------
    rsp : `~gammalib.GCTAResponse`
        Response function
    name : str
        Name of the response function
    color : str, optional
        Color for plot
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
    plt.xlabel('Energy (TeV)')
    plt.ylabel(r'Background rate (counts s$^{-1}$ MeV$^{-1}$ sr$^{-1}$)')

    # Set legend
    plt.legend(loc='lower left')

    # Return
    return


# ==================== #
# Show one sensitivity #
# ==================== #
def show_one_sensitivity(rsp, name, color='r', duration=180000.0, alpha=0.2,
                         sigma=5.0):
    """
    Show one sensitivity

    Parameters
    ----------
    rsp : `~gammalib.GCTAResponse`
        Response function
    name : str
        Name of the response function
    color : str, optional
        Color for plot
    duration : float, optional
        Duration of observation (s)
    alpha : float, optional
        Ratio of on-to-off exposure
    sigma : float, optional
        Required detection significance in Gaussian sigma
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

    # Initialise sensitivity vector
    flux = [0.0 for i in range(nbins)]

    # Loop over all energy bins
    for i in range(nbins):

        # Compute the energy width in MeV
        ewidth = (math.pow(10.0, logE[i]+0.1) - math.pow(10.0, logE[i]-0.1)) * \
                 1.0e6

        # Extract effective area in cm2 from Instrument Response Function
        aeff = rsp.aeff()(logE[i])

        # Extract 68% containment radius in radians from Instrument Response
        # Function.
        # A kluge is used since this information is not directly available,
        # but the "delta_max()" method returns 5 times the "sigma" of the
        # Gaussian PSF. Hence dividing by 5 gives the "sigma", and dividing
        # by the 68% containment radius to sigma conversion factor gives
        # r68.
        r68 = rsp.psf().delta_max(logE[i]) / 5.0 / r68_to_sigma

        # Compute the solid angle of the 68% containment radius
        solidangle = gammalib.twopi * (1.0 - math.cos(r68))

        # Compute the total number of background counts by multiplying the
        # background rate by the observation duration, the energy width and
        # the solid angle
        bgd_counts = rsp.background()(logE[i], 0.0, 0.0) * \
                     duration * ewidth * solidangle

        # The number of Off counts is the total number of background counts
        # divided by the "alpha" parameter
        Noff = bgd_counts/alpha

        # If there are Off counts then compute the corresponding On counts
        # for "sigma" detection significance
        if Noff > 0:

            # Compute On counts using Li & Ma
            Non = Non_lima(sigma, Noff)

            # Get the source counts by subtracting the  background counts
            # under the source. Make sure that the source counts are at least
            # 5% of the background counts (low-energy background systematics)
            # and that there are at least 10 source events (high-energy limit)
            src_counts = Non - bgd_counts
            if src_counts < 0.05*bgd_counts:
                src_counts = 0.05*bgd_counts
            if src_counts < 10:
                src_counts = 10.0

            # Compute flux under the assumption that the source spectrum is
            # a power law with index of -2.6. The flux is the source counts
            # divided by the duration and 68% of the effective area times
            # a conversion factor. The conversion factor is given by
            #
            #                 conv = 1.6021765 * E^2 / F(E)
            #
            # where E is the mean energy of the bin and F(E) is the flux in
            # the bin
            emin    = gammalib.GEnergy(math.pow(10.0, logE[i]-0.1), 'TeV')
            emax    = gammalib.GEnergy(math.pow(10.0, logE[i]+0.1), 'TeV')
            epivot  = gammalib.GEnergy(math.pow(10.0, logE[i]), 'TeV')
            plaw    = gammalib.GModelSpectralPlaw(1.0e-6, -2.6, epivot)
            conv    = TeV2erg*E[i]*E[i]/plaw.flux(emin, emax)
            flux[i] = conv * src_counts / (duration * aeff*0.68)

    # Plot data
    plt.loglog(E, flux, color+'-', label=name)
    plt.loglog(E, flux, color+'o')

    # Set axes
    plt.xlabel('Energy (TeV)')
    plt.ylabel(r'Sensitivity (erg cm$^{-2}$ s$^{-1}$)')

    # Set legend
    plt.legend(loc='upper right')

    # Return
    return


# ==================== #
# Show one sensitivity #
# ==================== #
def show_one_response(rspname, dbname, name, rootdir=None, color='r'):
    """
    Show one response.

    Parameters
    ----------
    rspname : str
        Response name
    dbname : str
        Database name
    name : str
        Name of the response function
    rootdir : str, optional
        Response root directory
    color : str, optional
        Color for plot
    """
    # Set-up calibration database
    caldb = gammalib.GCaldb()
    if rootdir != None:
        caldb.rootdir(rootdir)
    if gammalib.dir_exists(dbname):
        caldb.rootdir(dbname)
    else:
        caldb.open('cta', dbname)

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


# ==================== #
# Show one sensitivity #
# ==================== #
def plot_response(plotfile):
    """
    Show response

    Parameters
    ----------
    plotfile : str
        Plot filename
    """
    # Create figures
    plt.figure(1)
    plt.title('Effective area')
    plt.figure(2)
    plt.title('Background rate')
    plt.figure(3)
    plt.title('Sensitivity')

    # Set response dictionary
    rsps = [{'dbname':  'prod2',
             'rspname': 'South_0.5h',
             'name':    'Prod2 (30 min)',
             'color':   'b'},
            {'dbname':  'prod2',
             'rspname': 'South_50h',
             'name':    'Prod2 (50 h)',
             'color':   'r'}]

    # Loop over all responses
    for rsp in rsps:

        # Get attributes
        dbname  = rsp['dbname']
        rspname = rsp['rspname']
        color   = rsp['color']
        name    = rsp['name']
        rootdir = rsp.setdefault('rootdir', None)

        # Show response
        show_one_response(rspname, dbname, name, rootdir=rootdir, color=color)

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ============= #
# Show response #
# ============= #
def show_response():
    """
    Show response
    """
    # Set usage string
    usage = 'show_response.py [-p plotfile]'

    # Set default options
    options = [{'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']

    # Show response
    plot_response(plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show response
    show_response()
