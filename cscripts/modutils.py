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


# ================================= #
# Returns attributes for ds9 string #
# ================================= #
def ds9attributes(model, free_color='green', fixed_color='magenta',
                  width=2, fontfamily='helvetica', fontsize=12,
                  fontweight='normal', fontslant='roman'):
    """
    Returns attributes for DS9 string

    Parameters
    ----------
    model : `~gammalib.GModel`
        Model
    free_color : str, optional
        Color for sources with free parameters (any DS9 color or hex code)
    fixed_color : str, optional
        Color for source without free parameters (any DS9 color or hex code)
    width : int, optional
        Line width for regions
    fontfamily : str, optional
        Font family for source labels (helvetica,times,courier)
    fontsize : int, optional
        Font size for source labels
    fontweight : str, optional
        Font weight for source labels (normal,bold)
    fontslant : str, optional
        Font slant for source labels (roman,italic)

    Returns
    -------
    attributes : str
        DS9 attribute string
    """
    # Determine region color. A model where all parameters are fixed will
    # get the color defined by "fixed_color", a model with at least one
    # free parameter will get the color defined by "free_color".
    color = fixed_color
    for par in model:
        if par.is_free():
            color = free_color
            break

    # Set DS9 attributes
    attributes = (' color=%s width=%d font="%s %d %s %s"' %
                  (color, width, fontfamily, fontsize, fontweight, fontslant))

    # Return attributes
    return attributes


# =========================== #
# Convert model to ds9 string #
# =========================== #
def model2ds9string(model, pnt_type='cross', pnt_mark_size=12,
                    show_labels=True, show_ext_type=True,
                    free_color='green', fixed_color='magenta',
                    width=2, fontfamily='helvetica', fontsize=12,
                    fontweight='normal', fontslant='roman'):
    """
    Converts model into a DS9 region string

    Parameters
    ----------
    model : `~gammalib.GModel`
        Model
    pnt_type : str, optional
        Marker type for point sources (circle,box,diamond,cross,x,arrow,boxcircle)
    pnt_mark_size : integer, optional
        Marker size for point sources
    show_labels : boolean, optional
        Add source labels?
    show_ext_type : boolean, optional
        Show type of extended model in source name?
    free_color : str, optional
        Color for sources with free parameters (any DS9 color or hex code)
    fixed_color : str, optional
        Color for source without free parameters (any DS9 color or hex code)
    width : int, optional
        Line width for regions
    fontfamily : str, optional
        Font family for source labels (helvetica,times,courier)
    fontsize : int, optional
        Font size for source labels
    fontweight : str, optional
        Font weight for source labels (normal,bold)
    fontslant : str, optional
        Font slant for source labels (roman,italic)

    Returns
    -------
    ds9string : str
        DS9 region string
    msg : str
        Error message
    """
    # Initialise DS9 region string and message string
    ds9string = ''
    msg       = ''

    # Retrieve model sky direction. The model is skipped in case that
    # it does not provide a sky direction.
    is_valid = True  
    try:
        modelpos = model.spatial().dir()
    except AttributeError:            
        msg = ('Skip model "%s" since it has no sky direction.\n' %
               model.name())
        is_valid = False

    # Continue only if sky direction was found
    if is_valid: 

        # Retrieve model type and name   
        modeltype = model.type()
        modelname = model.name()

        # Handle point source
        if modeltype == 'PointSource':

            # Append point with Right Ascension and Declination to the DS9
            # string. The type of the point is specified by the "pnt_type"
            # parameter, the size of the point by the "pnt_mark_size"
            # parameter
            ds9string += ('point(%.6f,%.6f) # point=%s %d' %
                          (modelpos.ra_deg(), modelpos.dec_deg(),
                           pnt_type, pnt_mark_size))

        # Handle extended sources    
        elif modeltype == "ExtendedSource":

            # Retrieve spatial model
            spatial   = model.spatial()
            classname = spatial.classname()

            # Handle radial sources
            if 'Radial' in classname:

                # Retrieve short name of model class (e.g. "Gauss"
                # or "Disk)
                shorttype = classname.split('Radial')[-1]

                # Handle Disk and Shell model
                if (classname == 'GModelSpatialRadialDisk' or
                    classname == 'GModelSpatialRadialShell'):
                    size = spatial.radius()

                # Handle Gauss Model
                elif classname == 'GModelSpatialRadialGauss':
                    size = spatial.sigma()

                # Skip if source is unknown
                else:
                    msg = ('Skip model "%s" since the radial model "%s" '
                           'is unknown.\n' % (model.name(), classname))
                    is_valid = False

                # Append circle to DS9 string
                if is_valid:
                    ds9string += ('circle(%.6f,%.6f,%.6f) #' %
                                  (modelpos.ra_deg(), modelpos.dec_deg(),
                                   size*3600.0))
        
            # Handle elliptical sources 
            elif 'Elliptical' in classname:

                # Retrieve short name and source size
                shorttype = classname.split('Elliptical')[-1]
                size1     = spatial.semimajor()
                size2     = spatial.semiminor()
                angle     = spatial.posangle()

                # Append ellipse to DS9 string
                ds9string += ('ellipse(%.6f,%.6f,%.6f,%.6f,%.6f) #' %
                              (modelpos.ra_deg(), modelpos.dec_deg(),
                               size1*3600.0, size2*3600.0, angle+90.0))
            
            # Skip if source is neither radial nor elliptical      
            else:
                msg = ('Skip model "%s" since the model "%s" is neither '
                       'a point source, a radial source, nor an elliptical '
                       'source.\n' % (model.name(), classname))
                is_valid = False
            
            # Add short model type to modelname
            if show_ext_type:
                modelname +=' ('+shorttype+')'

        # Add DS9 attributes
        if is_valid:
            ds9string += ds9attributes(model, free_color=free_color,
                                              fixed_color=fixed_color,
                                              width=width,
                                              fontfamily=fontfamily,
                                              fontsize=fontsize,
                                              fontweight=fontweight,
                                              fontslant=fontslant)
            if show_labels:
                ds9string += ' text={'+modelname+'}'
        else:
            ds9string = ''
    
    # Return ds9 and message strings
    return ds9string, msg


# ========================= #
# Save models into ds9 file #
# ========================= #
def models2ds9file(models, filename,
                   pnt_type='cross', pnt_mark_size=12,
                   show_labels=True, show_ext_type=True,
                   free_color='green', fixed_color='magenta',
                   width=2, fontfamily='helvetica', fontsize=12,
                   fontweight='normal', fontslant='roman'):
    """
    Save models in a DS9 region file

    Parameters
    ----------
    model : `~gammalib.GModel`
        Model
    filename : str
        DS9 region file name
    pnt_type : str, optional
        Marker type for point sources (circle,box,diamond,cross,x,arrow,boxcircle)
    pnt_mark_size : integer, optional
        Marker size for point sources
    show_labels : boolean, optional
        Add source labels?
    show_ext_type : boolean, optional
        Show type of extended model in source name?
    free_color : str, optional
        Color for sources with free parameters (any DS9 color or hex code)
    fixed_color : str, optional
        Color for source without free parameters (any DS9 color or hex code)
    width : int, optional
        Line width for regions
    fontfamily : str, optional
        Font family for source labels (helvetica,times,courier)
    fontsize : int, optional
        Font size for source labels
    fontweight : str, optional
        Font weight for source labels (normal,bold)
    fontslant : str, optional
        Font slant for source labels (roman,italic)

    Returns
    -------
    errors : str
        Error message
    """
    # Initialise error string
    errors = ''
    
    # Open file   
    f = open(filename, 'w')
    
    # Write coordinate system
    f.write('fk5\n')
     
    # Loop over models
    for model in models:
        
        # Continue only if point source or extended source model
        if (model.type() == 'PointSource' or
            model.type() == 'ExtendedSource'):
            line, msg = model2ds9string(model,
                                        pnt_type=pnt_type,
                                        pnt_mark_size=pnt_mark_size,
                                        show_labels=show_labels,
                                        free_color=free_color,
                                        fixed_color=fixed_color,
                                        width=width,
                                        fontfamily=fontfamily,
                                        fontsize=fontsize,
                                        fontweight=fontweight,
                                        fontslant=fontslant)
            if len(line):
                f.write(line+'\n')
            elif len(msg):
                errors += msg+'\n'
        
        # Logging for diffuse components    
        elif model.type() == 'DiffuseSource':
            errors += 'Skipping diffuse model "'+model.name()+'".\n'
        
        # Logging for background components 
        else:
            errors += 'Skipping background model "'+model.name()+'".\n'

    # Close file
    f.close()

    # Return error message string
    return errors
