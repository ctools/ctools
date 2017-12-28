#! /usr/bin/env python
# ==========================================================================
# Display residuals generated by csresspec
#
# Copyright (C) 2017- Luigi Tibaldo
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
    from matplotlib import gridspec

    plt.figure()
    plt.close()
except (ImportError, RuntimeError):
    print('This script needs the "matplotlib" module')
    sys.exit()


# ================================= #
# Fill counts/model/residual arrays #
# ================================= #
def fill_cmr(row, counts, model, resid, e_counts, e_resid, c_counts, c_model,
             c_resid, algorithm):
    """
    Helper function that fills in the counts, model, residuals lists and
    calculate relative errors for given table row
    :param row: int row number
    :param counts: list
    :param model: list
    :param resid: list
    :param e_counts: list
    :param e_resid: list
    :param c_counts: `~gammalib.GFits.table.column'
    :param c_model: `~gammalib.GFits.table.column'
    :param c_resid: `~gammalib.GFits.table.column'
    :param algorithm: string
    :return: counts, model, resid, e_counts, e_resid
    """
    counts.append(c_counts.real(row))
    model.append(c_model.real(row))
    resid.append(c_resid.real(row))

    # calculate count error
    err = math.sqrt(c_counts.real(row))
    if err == 1:
        # this prevents visualization problem in matplotlib
        err = 0.99
    else:
        pass
    e_counts.append(err)

    # calculate residual error
    if algorithm == 'SUB':
        e_resid.append(err)
    elif algorithm == 'SUBDIV':
        e_resid.append(err / c_model.real(row))
    elif algorithm == 'SUBDIVSQRT':
        e_resid.append(err / math.sqrt(c_model.real(row)))
    elif algorithm == 'SIGNIFICANCE':
        e_resid.append(1.)

    return counts, model, resid, e_counts, e_resid


# ============= #
# Plot spectrum #
# ============= #
def plot_residuals(filename, plotfile, hdu):
    """
    Plot spectrum

    Parameters
    ----------
    filename : str
        Name of spectrum FITS file
    plotfile : str
        Plot file name
    hdu      : int
        Number of observation to plot
    """
    # Read spectrum file    
    fits = gammalib.GFits(filename)
    table = fits.table(1 + hdu)
    c_emin = table['Emin']
    c_emax = table['Emax']
    c_counts = table['Counts']
    c_model = table['Model']
    c_resid = table['Residuals']
    try:
        c_counts_off = table['Counts_Off']
        c_model_off = table['Model_Off']
        c_resid_off = table['Residuals_Off']
        is_onoff = True
    except:
        is_onoff = False

    # Initialise arrays to be filled
    ebounds = []
    ed_engs = []
    eu_engs = []
    em_engs = []
    counts = []
    e_counts = []
    model = []
    resid = []
    e_resid = []
    counts_off = []
    e_counts_off = []
    model_off = []
    resid_off = []
    e_resid_off = []

    # Residual algorithm
    algorithm = table.header()['ALGORITHM'].string()

    # Loop over rows of the file
    nrows = table.nrows()
    # add first energy boundary
    ebounds.append(c_emin.real(0))
    for row in range(nrows):
        # boundaries
        ebounds.append(c_emax.real(row))
        # geometrical mean energy and errors
        emean = math.sqrt(c_emin.real(row) * c_emax.real(row))
        em_engs.append(emean)
        ed_engs.append(emean - c_emin.real(row))
        eu_engs.append(c_emax.real(row) - emean)
        # counts, model, residuals
        counts, model, resid, e_counts, e_resid = fill_cmr(row, counts, model,
                                                           resid, e_counts,
                                                           e_resid, c_counts,
                                                           c_model, c_resid,
                                                           algorithm)
        # onoff
        if is_onoff:
            counts_off, model_off, resid_off, e_counts_off, e_resid_off = fill_cmr(
                row, counts_off,
                model_off, resid_off,
                e_counts_off,
                e_resid_off,
                c_counts_off,
                c_model_off, c_resid_off,
                algorithm)

    # Add model value to be compatible with plt.step
    model = [model[0]] + model
    if is_onoff:
        model_off = [model_off[0]] + model_off

    # Plot
    axarr = []
    if is_onoff:
        f = plt.figure(figsize=(8, 4))
        gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1])
        gs.update(hspace=0, left=0.09, right=0.97)
        axarr.append(f.add_subplot(gs[0]))
        axarr.append(f.add_subplot(gs[2]))
        axarr.append(f.add_subplot(gs[1], sharex=axarr[0]))
        axarr.append(f.add_subplot(gs[3], sharex=axarr[1]))
    else:
        f = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        gs.update(hspace=0)
        axarr.append(f.add_subplot(gs[0]))
        axarr.append(f.add_subplot(gs[1], sharex=axarr[0]))
    plt.setp(axarr[1].get_xticklabels(), visible=False)
    axarr[0].set_xscale('log')
    axarr[1].set_xscale('log')
    axarr[0].set_yscale('log')
    # Counts and model
    axarr[0].errorbar(em_engs, counts, yerr=e_counts, xerr=[ed_engs, eu_engs],
                      fmt='ko', capsize=0, linewidth=2, zorder=2, label='Data')
    axarr[0].step(ebounds, model, color='0.5', linewidth=2, zorder=1,
                  label='Model')
    axarr[0].set_ylabel('Counts')
    # Residuals
    axarr[1].errorbar(em_engs, resid, yerr=e_resid, xerr=[ed_engs, eu_engs],
                      fmt='ko', capsize=0, linewidth=2, zorder=2)
    axarr[1].axhline(0, color='0.5', linestyle='--')
    axarr[1].set_xlabel('Energy (TeV)')
    if algorithm == 'SUB':
        axarr[1].set_ylabel('Residuals (counts)')
    elif algorithm == 'SUBDIV':
        axarr[1].set_ylabel('Residuals (fraction)')
    elif algorithm == 'SUBDIVSQRT' or algorithm == 'SIGNIFICANCE':
        axarr[1].set_ylabel(r'Residuals ($\sigma$)')
    # On/Off
    if is_onoff:
        axarr[0].set_title('ON')
        axarr[2].set_title('OFF')
        plt.setp(axarr[3].get_xticklabels(), visible=False)
        axarr[2].set_xscale('log')
        axarr[3].set_xscale('log')
        axarr[2].set_yscale('log')
        # Counts and model
        axarr[2].errorbar(em_engs, counts_off, yerr=e_counts_off,
                          xerr=[ed_engs, eu_engs],
                          fmt='ko', capsize=0, linewidth=2, zorder=2,
                          label='Data')
        axarr[2].step(ebounds, model_off, color='0.5', linewidth=2, zorder=1,
                      label='Model')
        # Residuals
        axarr[3].errorbar(em_engs, resid_off, yerr=e_resid_off,
                          xerr=[ed_engs, eu_engs],
                          fmt='ko', capsize=0, linewidth=2, zorder=2)
        axarr[3].axhline(0, color='0.5', linestyle='--')
        axarr[3].set_xlabel('Energy (TeV)')

    # Add spectra of individual components
    skiplist = ['Counts', 'Model', 'Residuals',
                'Counts_Off', 'Model_Off', 'Residuals_Off',
                'Emin', 'Emax']
    for s in range(table.ncols()):
        if table[s].name() in skiplist:
            pass
        else:
            component = []
            for row in range(nrows):
                component.append(table[s].real(row))
            component = [component[0]] + component
            axarr[0].step(ebounds, component, zorder=0,
                          label=table[s].name())

    # Add legend
    axarr[0].legend(loc='best')

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    # Return
    return


# =============== #
# Show residuals  #
# =============== #
def show_residuals():
    """
    Show residuals
    """
    # Set usage string
    usage = 'show_spectrum.py [-p plotfile] [-h hdu] [file]'

    # Set default options
    options = [{'option': '-p', 'value': ''},
               {'option': '-h', 'value': 0}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']
    hdu = options[1]['value']

    # Show residuals
    plot_residuals(args[0], plotfile, hdu)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    # Show spectrum
    show_residuals()
