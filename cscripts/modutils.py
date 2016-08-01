# ==========================================================================
# Utility functions for model handling
#
# Copyright (C) 2016 Juergen Knoedlseder
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


# ============================================ #
# Return model for TS fitting of a test source #
# ============================================ #
def test_source(models, srcname, ra=None, dec=None, fitspat=False,
                fitspec=False):
    """
    Return model for TS fitting of a test source

    Parameters
    ----------
    models : `~gammalib.GModels`
        Input model container
    srcname : str
        Test source name
    ra : float, optional
        Right Ascension of test source (deg)
    dec : float, optional
        Declination of test source (deg)
    fitspat : bool, optional
        Fit spatial parameters of all models?
    fitspec : bool, optional
        Fit spectral parameters of all models?

    Returns
    -------
    model : `~gammalib.GModels`
        Model container for TS fitting
    """
    # Create a clone of the input model
    outmodels = models.copy()

    # Disable TS computation for all model components (will enable the
    # test source later)
    for model in outmodels:
        model.tscalc(False)

    # Get source model and enable TS computation
    model = outmodels[srcname]
    model.tscalc(True)

    # If source model has no "Prefactor" parameter then raise an exception
    if not model.has_par('Prefactor'):
        msg = ('Model "%s" has no parameter "Prefactor". Only spectral '
               'models with a "Prefactor" parameter are supported.' %
               srcname)
        raise RuntimeError(msg)

    # Set position of test source
    if ra != None and dec != None:
        if model.has_par('RA') and model.has_par('DEC'):
            model['RA'].value(ra)
            model['DEC'].value(dec)

    # Set possible spatial and spectral parameters
    spatial  = ['RA', 'DEC', 'Sigma', 'Radius', 'Width', 'PA',
                'MinorRadius', 'MajorRadius']
    spectral = ['Index', 'Index1', 'Index2', 'BreakEnergy', 'CutoffEnergy',
                'InverseCutoffEnergy']

    # Fit or fix spatial parameters
    for par in spatial:
        if model.has_par(par):
            if fitspat:
                model[par].free()
            else:
                model[par].fix()

    # Fit or fix spectral parameters
    for par in spectral:
        if model.has_par(par):
            if fitspec:
                model[par].free()
            else:
                model[par].fix()

    # Return model container
    return outmodels

