#! /usr/bin/env python
# ==========================================================================
# Display Instrument Response Function
#
# Copyright (C) 2017 Juergen Knoedlseder
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


# ========= #
# Plot Aeff #
# ========= #
def plot_aeff(sub, aeff, emin=None, emax=None, tmin=None, tmax=None,
              nengs=100, nthetas=100):
    """
    Plot effective area

    Parameters
    ----------
    sub : figure
        Subplot
    aeff : `~gammalib.GCTAAeff2d`
        Instrument Response Function
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    tmin : float, optional
        Minimum offset angle (deg)
    tmax : float, optional
        Maximum offset angle (deg)
    nengs : int, optional
        Number of energies
    nthetas : int, optional
        Number of offset angles
    """
    # Determine energy range
    ieng = aeff.table().axis('ENERG')
    neng = aeff.table().axis_bins(ieng)
    if emin == None:
        emin = aeff.table().axis_lo(ieng, 0)
    if emax == None:
        emax = aeff.table().axis_hi(ieng, neng-1)

    # Determine offset angle range
    itheta = aeff.table().axis('THETA')
    ntheta = aeff.table().axis_bins(itheta)
    if tmin == None:
        tmin = aeff.table().axis_lo(itheta, 0)
    if tmax == None:
        tmax = aeff.table().axis_hi(itheta, ntheta-1)

    # Use log energies
    emin = math.log10(emin)
    emax = math.log10(emax)

    # Set axes
    denergy     = (emax - emin)/(nengs-1)
    dtheta      = (tmax - tmin)/(nthetas-1)
    logenergies = [emin+i*denergy for i in range(nengs)]
    thetas      = [tmax-i*dtheta  for i in range(nthetas)]

    # Initialise image
    image = []

    # Loop over offset angles
    for theta in thetas:

        # Initialise row
        row = []

        # Loop over energies
        for logenergy in logenergies:

            # Get effective area value
            value = aeff(logenergy, theta*gammalib.deg2rad)

            # Append value
            row.append(value)

        # Append row
        image.append(row)

    # Plot image
    sub.imshow(image, extent=[emin,emax,tmin,tmax], aspect=0.5)

    # Show boundary contours
    contours = sub.contour(logenergies, thetas, image, [0.0], colors=('white'))
    sub.clabel(contours, inline=1, fontsize=8)

    # Plot title and axis
    sub.set_title('Effective area')
    sub.set_xlabel('log10(E/TeV)')
    sub.set_ylabel('Offset angle (deg)')

    # Return
    return


# ======== #
# Plot PSF #
# ======== #
def plot_psf(sub, psf, emin=None, emax=None, tmin=None, tmax=None,
             nengs=100, nthetas=100):
    """
    Plot Point Spread Function

    Parameters
    ----------
    sub : figure
        Subplot
    psf : `~gammalib.GCTAPsf2D`
        Instrument Response Function
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    tmin : float, optional
        Minimum offset angle (deg)
    tmax : float, optional
        Maximum offset angle (deg)
    nengs : int, optional
        Number of energies
    nthetas : int, optional
        Number of offset angles
    """
    # Determine energy range
    ieng = psf.table().axis('ENERG')
    neng = psf.table().axis_bins(ieng)
    if emin == None:
        emin = psf.table().axis_lo(ieng, 0)
    if emax == None:
        emax = psf.table().axis_hi(ieng, neng-1)

    # Determine offset angle range
    itheta = psf.table().axis('THETA')
    ntheta = psf.table().axis_bins(itheta)
    if tmin == None:
        tmin = psf.table().axis_lo(itheta, 0)
    if tmax == None:
        tmax = psf.table().axis_hi(itheta, ntheta-1)

    # Use log energies
    emin = math.log10(emin)
    emax = math.log10(emax)

    # Set axes
    denergy     = (emax - emin)/(nengs-1)
    dtheta      = (tmax - tmin)/(nthetas-1)
    logenergies = [emin+i*denergy for i in range(nengs)]
    thetas      = [tmax-i*dtheta  for i in range(nthetas)]

    # Initialise image
    image = []

    # Loop over offset angles
    for theta in thetas:

        # Initialise row
        row = []

        # Loop over energies
        for logenergy in logenergies:

            # Get containment radius value
            value = psf.containment_radius(0.68, logenergy, theta*gammalib.deg2rad) * \
                    gammalib.rad2deg

            # Append value
            row.append(value)

        # Append row
        image.append(row)

    # Plot image
    sub.imshow(image, extent=[emin,emax,tmin,tmax], aspect=0.5, vmin=0.0, vmax=0.5)

    # Show boundary contours
    contours = sub.contour(logenergies, thetas, image, [0.0], colors=('white'))
    sub.clabel(contours, inline=1, fontsize=8)

    # Plot title and axis
    if psf.classname() == 'GCTAPsfKing':
        sub.set_title('King function PSF 68% containment radius')
    else:
        sub.set_title('Gaussian PSF 68% containment radius')
    sub.set_xlabel('log10(E/TeV)')
    sub.set_ylabel('Offset angle (deg)')

    # Return
    return


# ====================== #
# Plot energy dispersion #
# ====================== #
def plot_edisp(sub, edisp, emin=None, emax=None, tmin=None, tmax=None,
               nengs=100, nthetas=100):
    """
    Plot Background template

    Parameters
    ----------
    sub : figure
        Subplot
    edisp : `~gammalib.GCTAEdisp2D`
        Instrument Response Function
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    tmin : float, optional
        Minimum offset angle (deg)
    tmax : float, optional
        Maximum offset angle (deg)
    nengs : int, optional
        Number of energies
    nthetas : int, optional
        Number of offset angles
    """
    # Determine energy range
    ieng = edisp.table().axis('ETRUE')
    neng = edisp.table().axis_bins(ieng)
    if emin == None:
        emin = edisp.table().axis_lo(ieng, 0)
    if emax == None:
        emax = edisp.table().axis_hi(ieng, neng-1)

    # Determine offset angle range
    itheta = edisp.table().axis('THETA')
    ntheta = edisp.table().axis_bins(itheta)
    if tmin == None:
        tmin = edisp.table().axis_lo(itheta, 0)
    if tmax == None:
        tmax = edisp.table().axis_hi(itheta, ntheta-1)

    # Use log energies
    emin = math.log10(emin)
    emax = math.log10(emax)

    # Set axes
    denergy     = (emax - emin)/(nengs-1)
    dtheta      = (tmax - tmin)/(nthetas-1)
    logenergies = [emin+i*denergy for i in range(nengs)]
    thetas      = [tmax-i*dtheta  for i in range(nthetas)]

    # Initialise image
    image = []

    # Loop over offset angles
    for theta in thetas:

        # Initialise row
        row = []

        # Compute detx and dety
        detx = theta*gammalib.deg2rad
        dety = 0.0

        # Loop over energies
        for logenergy in logenergies:

            # Get no migration energy dispersion value
            value = edisp(logenergy, logenergy, theta*gammalib.deg2rad)

            # Append value
            row.append(value)

        # Append row
        image.append(row)

    # Plot image
    sub.imshow(image, extent=[emin,emax,tmin,tmax], aspect=0.5)

    # Show boundary contours
    contours = sub.contour(logenergies, thetas, image, [0.0], colors=('white'))
    sub.clabel(contours, inline=1, fontsize=8)

    # Plot title and axis
    sub.set_title('No migration energy dispersion')
    sub.set_xlabel('log10(E/TeV)')
    sub.set_ylabel('Offset angle (deg)')

    # Return
    return


# =============== #
# Plot Background #
# =============== #
def plot_bkg(sub, bkg, emin=None, emax=None, tmin=None, tmax=None,
             nengs=100, nthetas=100):
    """
    Plot Background template

    Parameters
    ----------
    sub : figure
        Subplot
    bkg : `~gammalib.GCTABackground3D`
        Instrument Response Function
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    tmin : float, optional
        Minimum offset angle (deg)
    tmax : float, optional
        Maximum offset angle (deg)
    nengs : int, optional
        Number of energies
    nthetas : int, optional
        Number of offset angles
    """
    # Determine energy range
    ieng = bkg.table().axis('ENERG')
    neng = bkg.table().axis_bins(ieng)
    if emin == None:
        emin = bkg.table().axis_lo(ieng, 0)
    if emax == None:
        emax = bkg.table().axis_hi(ieng, neng-1)

    # Determine offset angle range
    if tmin == None:
        tmin = 0.0
    if tmax == None:
        tmax = 6.0

    # Use log energies
    emin = math.log10(emin)
    emax = math.log10(emax)

    # Set axes
    denergy     = (emax - emin)/(nengs-1)
    dtheta      = (tmax - tmin)/(nthetas-1)
    logenergies = [emin+i*denergy for i in range(nengs)]
    thetas      = [tmax-i*dtheta  for i in range(nthetas)]

    # Initialise image
    image = []

    # Loop over offset angles
    for theta in thetas:

        # Initialise row
        row = []

        # Compute detx and dety
        detx = theta*gammalib.deg2rad
        dety = 0.0

        # Loop over energies
        for logenergy in logenergies:

            # Get containment radius value
            value = bkg(logenergy, detx, dety)

            # Append value
            row.append(value)

        # Append row
        image.append(row)

    # Plot image
    sub.imshow(image, extent=[emin,emax,tmin,tmax], aspect=0.5)

    # Show boundary contours
    contours = sub.contour(logenergies, thetas, image, [0.0], colors=('white'))
    sub.clabel(contours, inline=1, fontsize=8)

    # Plot title and axis
    sub.set_title('Background acceptance')
    sub.set_xlabel('log10(E/TeV)')
    sub.set_ylabel('Offset angle (deg)')

    # Return
    return


# ======== #
# Plot IRF #
# ======== #
def plot_irf(irf, emin, emax, tmin, tmax, plotfile):
    """
    Plot Instrument Response Function

    Parameters
    ----------
    irf : `~gammalib.GCTAResponseIrf`
        Instrument Response Function
    emin : float
        Minimum energy (TeV)
    emax : float
        Maximum energy (TeV)
    tmin : float
        Minimum offset angle (deg)
    tmax : float
        Maximum offset angle (deg)
    plotfile : str
        Plot filename
    """
    # Build selection string
    selection  = ''
    eselection = ''
    tselection = ''
    if emin != None and emax != None:
        eselection += '%.3f-%.1f TeV' % (emin, emax)
    elif emin != None:
        eselection += ' >%.3f TeV' % (emin)
    elif emax != None:
        eselection += ' <%.1f TeV' % (emax)
    if tmin != None and tmax != None:
        tselection += '%.1f-%.1f deg' % (tmin, tmax)
    elif tmin != None:
        tselection += ' >%.1f deg' % (tmin)
    elif tmax != None:
        tselection += ' <%.1f deg' % (tmax)
    if len(eselection) > 0 and len(tselection) > 0:
        selection = ' (%s, %s)' % (eselection, tselection)
    elif len(eselection) > 0:
        selection = ' (%s)' % (eselection)
    elif len(tselection) > 0:
        selection = ' (%s)' % (tselection)

    # Build title
    mission    = irf.caldb().mission()
    instrument = irf.caldb().instrument()
    response   = irf.rspname()
    title      = '%s "%s" Instrument Response Function "%s"%s' % \
                 (gammalib.toupper(mission), instrument, response, selection)

    # Create figure
    fig = plt.figure(figsize=(12,8))

    # Add title
    fig.suptitle(title, fontsize=16)

    # Plot Aeff
    ax1 = fig.add_subplot(221)
    plot_aeff(ax1, irf.aeff(), emin=emin, emax=emax, tmin=tmin, tmax=tmax)

    # Plot Psf
    ax2 = fig.add_subplot(222)
    plot_psf(ax2, irf.psf(), emin=emin, emax=emax, tmin=tmin, tmax=tmax)

    # Plot Edisp
    ax3 = fig.add_subplot(223)
    plot_edisp(ax3, irf.edisp(), emin=emin, emax=emax, tmin=tmin, tmax=tmax)

    # Plot Background
    ax4 = fig.add_subplot(224)
    plot_bkg(ax4, irf.background(), emin=emin, emax=emax, tmin=tmin, tmax=tmax)

    # Show plots or save it into file
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# ======== #
# Show IRF #
# ======== #
def show_irf():
    """
    Show Instrument Response Function
    """
    # Set usage string
    usage = 'show_irf.py [-p plotfile] caldb irf'

    # Set default options
    options = [{'option': '-p',    'value': ''},
               {'option': '-emin', 'value': None},
               {'option': '-emax', 'value': None},
               {'option': '-tmin', 'value': None},
               {'option': '-tmax', 'value': None}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']
    emin     = options[1]['value']
    emax     = options[2]['value']
    tmin     = options[3]['value']
    tmax     = options[4]['value']

    # Convert limits to float
    if emin != None:
        emin = float(emin)
    if emax != None:
        emax = float(emax)
    if tmin != None:
        tmin = float(tmin)
    if tmax != None:
        tmax = float(tmax)

    # Get IRF
    caldb = gammalib.GCaldb('cta', args[0])
    irf   = gammalib.GCTAResponseIrf(args[1], caldb)
    
    # Plot IRF
    plot_irf(irf, emin, emax, tmin, tmax, plotfile)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show IRF
    show_irf()
