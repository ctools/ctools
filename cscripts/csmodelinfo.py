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
import gammalib
import ctools
import sys


# ================= #
# csmodelinfo class #
# ================= #
class csmodelinfo(ctools.cscript):
    """
    Shows the content of a model container.

    The script is useful to show, e.g., the number of sources in a model
    container, the background models, free and fixed parameters etc.
    The script can also dump the models into a region file to for display
    with DS9.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name    = "csmodelinfo"
        self._version = "1.1.0"

        # Initialise class members
        self._models        = gammalib.GModels()
        self._ds9file       = gammalib.GFilename("NONE")
        self._pnt_type      = ""
        self._pnt_mark_size = 12
        self._show_labels   = True
        self._width         = 2
        self._fontfamily    = "helvetica"
        self._fontsize      = 12
        self._fontweight    = "normal"
        self._fontslant     = "roman"
        self._show_ext_type = True
        self._free_color    = "green"
        self._fixed_color   = "magenta"

        # Initialise application by calling the appropriate class
        # constructor.
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get models
        self._models = gammalib.GModels(self["inmodel"].filename())
        
        # Get hidden parameters for region file
        self._pnt_type      = self["pnt_type"].string()
        self._pnt_mark_size = self["pnt_mark_size"].integer()
        self._show_labels   = self["show_labels"].boolean()
        self._width         = self["width"].integer()
        self._fontfamily    = self["fontfamily"].string()
        self._fontsize      = self["fontsize"].integer()
        self._fontweight    = self["fontweight"].string()
        self._fontslant     = self["fontslant"].string()
        self._show_ext_type = self["show_ext_type"].boolean()
        self._free_color    = self["free_color"].string()
        self._fixed_color   = self["fixed_color"].string()

        # Read ahead DS9 filename
        if self._read_ahead():
            self._ds9file = self["ds9file"].filename()

        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")
        
        # Return
        return

    def _model2ds9string(self, model):
        """
        Converts model into a DS9 region string.

        Args:
            model: Model.

        Returns:
            DS9 region string.
        """
        # Initialise DS9 region string
        ds9string = ""

        # Determine region color. The color will change between a model
        # where all parameters are fixed and a model with at least one
        # free parameter.
        color = self._fixed_color
        for par in model:
            if par.is_free():
                color = self._free_color
                break

        # Retrieve model sky direction. The model is skipped in case that
        # it does not provide a sky direction.
        is_valid = True  
        try:
            modelpos = model.spatial().dir()
        except AttributeError:            
            self._log("Skip source \""+model.name()+"\" since no "
                      "sky direction is associated with that model.\n")
            is_valid = False

        # Continue only if sky direction was found
        if is_valid: 

            # Retrieve model type and name   
            modeltype = model.type()
            modelname = model.name()

            # Handle point source
            if modeltype == "PointSource":

                # Append point to DS9 string
                ds9string += "point("
                ds9string += str(modelpos.ra_deg()) 
                ds9string += ","
                ds9string += str(modelpos.dec_deg())
                ds9string += ") # "
                ds9string += "point="+self._pnt_type+" "
                ds9string += str(self._pnt_mark_size)+" "

            # Handle extended sources    
            elif modeltype == "ExtendedSource":

                # Retrieve spatial model
                spatial   = model.spatial()
                classname = spatial.classname()

                # Handle radial sources
                if "Radial" in classname:  

                    # Retrieve short name of model class (e.g. "Gauss"
                    # or "Disk)
                    shorttype = classname.split("Radial")[-1]

                    # Handle Disk and Shell model
                    if (classname == "GModelSpatialRadialDisk" or
                        classname == "GModelSpatialRadialShell"):
                        size = spatial.radius()

                    # Handle Gauss Model
                    elif classname == "GModelSpatialRadialGauss":
                        size = spatial.sigma()

                    # Skip if source is unknown
                    else:
                        self._log("Skip source \""+model.name()+"\""
                                  " since the radial model \""+classname+
                                  "\" is unkown.\n")
                        is_valid = False

                    # Append circle to DS9 string
                    if is_valid:
                        ds9string += "circle("
                        ds9string += str(modelpos.ra_deg()) 
                        ds9string += ","
                        ds9string += str(modelpos.dec_deg())
                        ds9string += ","
                        ds9string += str(size*3600.)+"\") # "
                
                # Handle elliptical sources 
                elif "Elliptical" in classname:

                    # Retrieve short name and source size
                    shorttype = classname.split("Elliptical")[-1]
                    size1     = spatial.semimajor()
                    size2     = spatial.semiminor()
                    angle     = spatial.posangle()

                    # Append ellipse to DS9 string
                    ds9string += "ellipse("
                    ds9string += str(modelpos.ra_deg()) 
                    ds9string += ","
                    ds9string += str(modelpos.dec_deg())
                    ds9string += ","
                    ds9string += str(size1*3600.)+"\""
                    ds9string += ","
                    ds9string += str(size2*3600.)+"\""
                    ds9string += ","
                    ds9string += str(angle+90.0)+") # "   
                
                # Skip if source is neither radial nor elliptical      
                else:
                    self._log("Skip source \""+model.name()+"\" "
                              "since the model \""+classname+"\" is "
                              "neither a point source, a radial source, "
                              "nor a elliptical source.\n")
                    is_valid = False
                
                # Add short model type to modelname
                if self._show_ext_type:
                    modelname +=" ("+shorttype+")"
    
            # Add DS9 attributes
            if is_valid:
                ds9string += "color="+color+" "
                ds9string += "width="+str(self._width)+" "
                ds9string += "font=\""+self._fontfamily+" "
                ds9string += str(self._fontsize)+" "
                ds9string += self._fontweight+" "
                ds9string += self._fontslant+"\""
                if self._show_labels:
                    ds9string += " text={"+modelname+"}"
        
        # Return string
        return ds9string 
        

    # Public methods
    def run(self):
        """
        Run the script.
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()
        
        # Optionally log models
        if self._logVerbose():
            self._log("\n")
            self._log.header1("Models")
            self._log("\n")
            self._log(str(self._models))
            self._log("\n")
        
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
            skymodel = (model.classname() == "GModelSky")
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
            self._log("\n")
            self._log.header1("Summary")
        
            # Log instruments
            self._log.header3("Instrument specific models")
            for inst, n_inst in instruments.items():
                if inst == "":
                    inst = "All"
                self._log.parformat(inst)
                self._log(str(n_inst))
                self._log("\n")
        
            # Log model types
            self._log.header3("Model types")
            for modeltype, n_types in types.items():
                self._log.parformat(modeltype)
                self._log(str(n_types))
                self._log("\n")
        
            # Log model parameter information
            self._log("\n")
            self._log.header1("Parameter information")
            self._log.parformat("All parameters")
            self._log(str(n_par_total))
            self._log("\n")
            self._log.parformat("Fixed parameters")
            self._log(str(n_par_fixed))
            self._log("\n")
            self._log.parformat("Free parameters (total)")
            self._log(str(n_par_free))
            self._log("\n")
            self._log.parformat("Free background parameters")
            self._log(str(n_par_free_bkg))
            self._log("\n")
            self._log.parformat("Free source parameters")
            self._log(str(n_par_free_src))
            self._log("\n")
            if n_par_free_spec > 0:
                self._log.parformat("Free spectral parameters")
                self._log(str(n_par_free_spec))
                self._log("\n")
            if n_par_free_spat > 0:
                self._log.parformat("Free spatial parameters")
                self._log(str(n_par_free_spat))
                self._log("\n")
            if n_par_free_temp > 0:
                self._log.parformat("Free temporal parameters")
                self._log(str(n_par_free_temp))
                self._log("\n")
            self._log.parformat("Parameters at limit")
            self._log(str(n_par_at_limit))
            self._log("\n")
        
            # Log parameters at limit (if any)
            if n_par_at_limit > 0:
                self._log.header3("Parameters at limit")
                for source, parameter in pars_at_limit.items():   
                    if len(parameter):  
                        for par in parameter:
                            self._log.parformat(source)
                            self._log(par)
                            self._log("\n")
        
            # Optionally log free parameters
            if self._logExplicit() and len(free_src_pars):
                self._log("\n")
                self._log.header2("Free source parameters")
                for source, parameter in free_src_pars.items():   
                    if len(parameter): 
                        self._log.header3(source)
                        for par in parameter:
                            self._log(str(par))
                            self._log("\n")
                        self._log("\n")            
        
            # Log TS values if available
            if len(ts):       
                self._log.header3("Test statistics")
                for source,tsvalue in ts.items():     
                    self._log.parformat(source)
                    self._log("%.2f"%tsvalue)
                    self._log("\n")     
                
        # Return
        return

    def save(self):
        """ 
        Save models to ds9 region file if required
        """
        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Save models in DS9 file")

        # Get output filename in case it was not read ahead
        self._ds9file = self["ds9file"].filename()

        # Check if DS9 file is valid
        if self._ds9file != "NONE":      
            
            # Log filename
            if self._logTerse():
                self._log(gammalib.parformat("DS9 filename"))
                self._log(self._ds9file.url())
                self._log("\n")
            
            # Open file   
            f = open(self._ds9file.url(),"w")
            
            # Write coordinate system
            f.write("fk5\n")
             
            # Loop over models
            for model in self._models:
                
                # Continue only if point source or extended source model
                if (model.type() == "PointSource" or
                    model.type() == "ExtendedSource"):   
                    line = self._model2ds9string(model)
                    if len(line):
                        f.write(line+"\n") 
                
                # Logging for diffuse components    
                elif model.type() == "DiffuseSource":
                    self._log("Skipping diffuse model \""+model.name()+
                              "\"\n")
                    
                # Logging for background components 
                else:
                    if self._logExplicit():
                        self._log("Skipping background model \""+
                                  model.name()+"\"\n")   
                                      
            # Close file
            f.close()
        
        # Return
        return 

    def execute(self):
        """
        Execute the script.
        """
        # Open logfile
        self._logFileOpen()

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
