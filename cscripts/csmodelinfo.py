#!/usr/bin/env python
# ==========================================================================
# Shows the content of a model container
#
# Copyright (C) 2015-2016 Michael Mayer
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
import gammalib
import ctools
from cscripts import modutils


# ================= #
# csmodelinfo class #
# ================= #
class csmodelinfo(ctools.cscript):
    """
    Shows the content of a model container

    The script is useful to show, e.g., the number of sources in a model
    container, the background models, free and fixed parameters etc.
    The script can also dump the models into a region file to for display
    with DS9.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name
        self._name    = 'csmodelinfo'
        self._version = '1.2.0'

        # Initialise class members
        self._models = gammalib.GModels()

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation
        """
        # Get models
        self._models = gammalib.GModels(self['inmodel'].filename())
        
        # Query hidden parameters for region file
        self['pnt_type'].string()
        self['pnt_mark_size'].integer()
        self['show_labels'].boolean()
        self['width'].integer()
        self['fontfamily'].string()
        self['fontsize'].integer()
        self['fontweight'].string()
        self['fontslant'].string()
        self['show_ext_type'].boolean()
        self['free_color'].string()
        self['fixed_color'].string()

        # Query ahead DS9 filename
        if self._read_ahead():
            self['outds9file'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)
        
        # Return
        return


    # Public methods
    def run(self):
        """
        Run the script
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()
        
        # Optionally log models
        if self._logVerbose():
            self._log('\n')
            self._log.header1('Models')
            self._log('\n')
            self._log(str(self._models))
            self._log('\n')
        
        # Initialise output values 
        types           = {}
        instruments     = {}
        ts              = {}
        pars_at_limit   = {}
        free_src_pars   = {}
        n_par_total     = 0
        n_par_free      = 0
        n_par_fixed     = 0
        n_par_at_limit  = 0
        n_par_free_bkg  = 0
        n_par_free_src  = 0
        n_par_free_spec = 0
        n_par_free_spat = 0
        n_par_free_temp = 0
        
        # Loop over models and initialise dictionaries with model names
        for model in self._models:
            types[model.type()] = 0 
            instruments[model.instruments()] = 0
            pars_at_limit[model.name()] = []

        # Loop over models again and fill values to be dumped
        for model in self._models:
            
            # Model types
            types[model.type()] += 1
            
            # Model instruments
            instruments[model.instruments()] += 1 
            
            # TS
            if model.has_ts():
                ts[model.name()] = model.ts()
                
            # Skymodel?
            skymodel = (model.classname() == 'GModelSky')
            if skymodel:
                free_src_pars[model.name()] = []                
            
            # Loop over model parameters
            for par in model:
                
                # Increment total parameters
                n_par_total += 1
                
                # Check if parameter is at limit
                if par.has_range():
                    if par.value() == par.min() or par.value() == par.max():
                        pars_at_limit[model.name()].append(par.name())   
                        n_par_at_limit += 1 
                        
                # Check if model parameter is free/fixed and
                # increment respective variable      
                if par.is_free():         
                    n_par_free += 1
                    if skymodel:
                        free_src_pars[model.name()].append(par)
                        n_par_free_src += 1
                    else:
                        n_par_free_bkg += 1
                else:
                    n_par_fixed += 1
            
            # If sky model, divide free parameters into spatial, spectral
            # and temporal
            if skymodel:
                for par in model.spectral():
                    if par.is_free():
                        n_par_free_spec += 1
                for par in model.spatial():
                    if par.is_free():
                        n_par_free_spat += 1
                for par in model.temporal():
                    if par.is_free():
                        n_par_free_temp += 1     

        # Optionally log model information
        if self._logTerse():
        
            # Log summary
            self._log_header1(gammalib.TERSE, 'Summary')
        
            # Log instruments
            self._log.header3('Instrument specific models')
            for inst, n_inst in instruments.items():
                if inst == '':
                    inst = 'All'
                self._log_value(gammalib.TERSE, inst, n_inst)
        
            # Log model types
            self._log_header3(gammalib.TERSE, 'Model types')
            for modeltype, n_types in types.items():
                self._log_value(gammalib.TERSE, modeltype, n_types)
        
            # Log model parameter information
            self._log_header1(gammalib.TERSE,'Parameter information')
            self._log_value(gammalib.TERSE,'All parameters', n_par_total)
            self._log_value(gammalib.TERSE,'Fixed parameters', n_par_fixed)
            self._log_value(gammalib.TERSE,'Free parameters (total)', n_par_free)
            self._log_value(gammalib.TERSE,'Free background parameters', n_par_free_bkg)
            self._log_value(gammalib.TERSE,'Free source parameters', n_par_free_src)
            if n_par_free_spec > 0:
                self._log_value(gammalib.TERSE,'Free spectral parameters', n_par_free_spec)
            if n_par_free_spat > 0:
                self._log_value(gammalib.TERSE,'Free spatial parameters', n_par_free_spat)
            if n_par_free_temp > 0:
                self._log_value(gammalib.TERSE,'Free temporal parameters', n_par_free_temp)
            self._log_value(gammalib.TERSE,'Parameters at limit', n_par_at_limit)
        
            # Log parameters at limit (if any)
            if n_par_at_limit > 0:
                self._log_header3(gammalib.TERSE, 'Parameters at limit')
                for source, parameter in pars_at_limit.items():   
                    if len(parameter):  
                        for par in parameter:
                            self._log_value(gammalib.TERSE, source, par)
        
            # Optionally log free parameters
            if len(free_src_pars):
                self._log_header2(gammalib.EXPLICIT, 'Free source parameters')
                for source, parameter in free_src_pars.items():   
                    if len(parameter): 
                        self._log.header3(source)
                        for par in parameter:
                            self._log_string(gammalib.EXPLICIT, str(par)+'\n')
        
            # Log TS values if available
            if len(ts):       
                self._log_header3(gammalib.TERSE, 'Test statistics')
                for source, tsvalue in ts.items():
                    self._log_value(gammalib.TERSE, source, tsvalue)
                
        # Return
        return

    def save(self):
        """ 
        Save models to ds9 region file
        """
        # Get output filename
        ds9file = self['outds9file'].filename()

        # Check if DS9 file is valid
        if ds9file != '' and ds9file != 'NONE':
            
            # Write header
            self._log_header1(gammalib.TERSE, 'Save models in DS9 file')
            
            # Log filename
            self._log_value(gammalib.NORMAL, 'DS9 filename', ds9file.url())

            # Get ds9 parameters
            pnt_type      = self['pnt_type'].string()
            pnt_mark_size = self['pnt_mark_size'].integer()
            show_labels   = self['show_labels'].boolean()
            free_color    = self['free_color'].string()
            fixed_color   = self['fixed_color'].string()
            width         = self['width'].integer()
            fontfamily    = self['fontfamily'].string()
            fontsize      = self['fontsize'].integer()
            fontweight    = self['fontweight'].string()
            fontslant     = self['fontslant'].string()

            # Save DS9 region file
            errors = modutils.models2ds9file(self._models, ds9file.url(),
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

            # Log errors
            if len(errors):
                self._log(errors)
    
        # Return
        return 

    def execute(self):
        """
        Execute the script
        """
        # Open logfile
        self.logFileOpen()

        # Read ahead output parameters
        self._read_ahead(True)

        # Run the script
        self.run()

        # Save ds9 file if required
        self.save()

        # Return
        return    
        

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csmodelinfo(sys.argv)
    
    # Execute application
    app.execute()
