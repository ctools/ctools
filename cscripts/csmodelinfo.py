#!/usr/bin/env python
# ==========================================================================
# Dump information about models into log file
#
# Copyright (C) 2015 Michael Mayer
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
    This class dumps information about a model container into a logfile or
    on screen. This might be helpful for quick access to the model container,
    showing, e.g., number of sources, background models, free and fixed
    parameters etc. The tool can also dump the models into a region file to
    display e.g. with ds9.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "csmodelinfo"
        self.version = "1.1.0"

        # Make sure that parfile exists
        file = self.parfile()

        # Initialise application
        if len(argv) == 0:
            ctools.cscript.__init__(self, self.name, self.version)
        elif len(argv) == 1:
            ctools.cscript.__init__(self, self.name, self.version, *argv)
        else:
            raise TypeError("Invalid number of arguments given.")

        # Set logger properties
        self.log_header()
        self.log.date(True)

        # Return
        return

    def __del__(self):
        """
        Destructor.
        """
        # Return
        return

    def parfile(self):
        """
        Check if parfile exists. If parfile does not exist then create a
        default parfile. This kluge avoids shipping the cscript with a
        parfile.
        """

        # Set parfile name
        parfile = self.name+".par"

        try:
            pars = gammalib.GApplicationPars(parfile)
        except:
            # Signal that parfile was not found
            print("Parfile \""+parfile+"\" not found. Create default parfile.")

            # Create default parfile
            pars = gammalib.GApplicationPars()
            pars.append(gammalib.GApplicationPar("inmodel","f","a","model.xml","","","Input model XML file"))
            pars.append(gammalib.GApplicationPar("ds9file","f","h","NONE","","","DS9 region file containing models"))
            pars.append(gammalib.GApplicationPar("pnt_type","s","h","cross","circle|box|diamond|cross|x|arrow|boxcircle","","Marker type for point sources"))
            pars.append(gammalib.GApplicationPar("pnt_mark_size","i","h","12","","","Marker size for point sources"))
            pars.append(gammalib.GApplicationPar("show_labels","b","h","yes","yes|no","","Add source labels?"))
            pars.append(gammalib.GApplicationPar("width","i","h","2","","","Line width for regions"))
            pars.append(gammalib.GApplicationPar("font","s","h","helvetica","helvetica|times|courier","","Font for source labels"))
            pars.append(gammalib.GApplicationPar("fontsize","i","h","12","","","Font size for source labels"))
            pars.append(gammalib.GApplicationPar("fonttype","s","h","normal","normal|bold","","Use normal or bold font?"))
            pars.append(gammalib.GApplicationPar("fonttype2","s","h","roman","roman|italic","","Use roman or italic font?"))
            pars.append(gammalib.GApplicationPar("show_ext_type","b","h","yes","yes|no","","Show type of extended model in source name"))
            pars.append(gammalib.GApplicationPar("free_color","s","h","green","","","Color for sources with free parameters (any ds9 color or hex code)"))
            pars.append(gammalib.GApplicationPar("fixed_color","s","h","magenta","","","Color for source without free parameters (any ds9 color or hex code)"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csmodelinfo.log","","","Log filename"))
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get models
        self.models = gammalib.GModels(self["inmodel"].filename())
        
        # Get ds9 filename
        self.ds9file = self["ds9file"].filename()
        
        # Get hidden parameters for region file
        self.pnt_type      = self["pnt_type"].string()
        self.pnt_mark_size = self["pnt_mark_size"].integer()
        self.show_labels   = self["show_labels"].boolean()
        self.width         = self["width"].integer()
        self.font          = self["font"].string()
        self.fontsize      = self["fontsize"].integer()
        self.fonttype      = self["fonttype"].string()
        self.fonttype2     = self["fonttype2"].string()
        self.show_ext_type = self["show_ext_type"].boolean()
        self.free_color    = self["free_color"].string()
        self.fixed_color   = self["fixed_color"].string()
        
        # Return
        return
        
    def run(self):
        """
        Run the script.
        """
        # Switch screen logging on in debug mode
        if self.logDebug():
            self.log.cout(True)

        # Get parameters
        self.get_parameters()
        
        # Write input parameters into logger
        if self.logTerse():
            self.log_parameters()
            self.log("\n")
            
        if self.logVerbose():
            self.log("\n")
            self.log.header1("Models in detail")
            self.log("\n")
            self.log(str(self.models))
            self.log("\n")
        
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
        for model in self.models:
            types[model.type()] = 0 
            instruments[model.instruments()] = 0
            pars_at_limit[model.name()] = []

        # Loop over models again and fill values to be dumped
        for model in self.models:
            
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
            
            # If sky model, divide free parameters into spatial, spectral and temporal
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
        
        # Log summary
        self.log("\n")
        self.log.header1("Summary")
        self.log("\n")
        
        # Log instruments
        self.log.header3("Instrument specific models")
        for inst, n_inst in instruments.iteritems():
            if inst == "":
                inst = "All"
            self.log.parformat(inst)
            self.log(str(n_inst))
            self.log("\n")
        
        # Log model types
        self.log("\n")
        self.log.header3("Model types")
        for modeltype, n_types in types.iteritems():
            self.log.parformat(modeltype)
            self.log(str(n_types))
            self.log("\n")
        
        # Log parameter information
        self.log("\n")
        self.log.header2("Parameter information")
        self.log.parformat("All parameters")
        self.log(str(n_par_total))
        self.log("\n")
        self.log.parformat("Fixed parameters")
        self.log(str(n_par_fixed))
        self.log("\n")
        self.log.parformat("Free parameters (total)")
        self.log(str(n_par_free))
        self.log("\n")
        self.log.parformat("Free background parameters")
        self.log(str(n_par_free_bkg))
        self.log("\n")
        self.log.parformat("Free source parameters")
        self.log(str(n_par_free_src))
        self.log("\n")
        if n_par_free_spec > 0:
            self.log.parformat("Free spectral parameters")
            self.log(str(n_par_free_spec))
            self.log("\n")
        if n_par_free_spat > 0:
            self.log.parformat("Free spatial parameters")
            self.log(str(n_par_free_spat))
            self.log("\n")
        if n_par_free_temp > 0:
            self.log.parformat("Free temporal parameters")
            self.log(str(n_par_free_temp))
            self.log("\n")
        self.log.parformat("Parameters at limit")
        self.log(str(n_par_at_limit))
        self.log("\n")
        self.log("\n")
        
        # Log parameters at limit (if any)
        if n_par_at_limit > 0:
            self.log.header3("Parameters at limit")
            for source, parameter in pars_at_limit.iteritems():   
                if len(parameter):  
                    for par in parameter:
                        self.log.parformat(source)
                        self.log(par)
                        self.log("\n")
        
        # Log free parameters if requested
        if self.logExplicit() and len(free_src_pars):
            self.log("\n")
            self.log.header2("Free source parameters")
            for source, parameter in free_src_pars.iteritems():   
                if len(parameter): 
                    self.log.header3(source)
                    for par in parameter:
                        self.log(str(par))
                        self.log("\n")
                    self.log("\n")            
        
        # Log TS values if available
        if len(ts):       
            self.log.header3("Test statistics")
            for source,tsvalue in ts.iteritems():     
                self.log.parformat(source)
                self.log("%.2f"%tsvalue)
                self.log("\n")     
                
        # Return
        return

    def save(self):
        """ 
        Save models to ds9 region file if required
        """
       
        # Check if ds9 file is valid
        if not self.ds9file == "NONE":      
            
            self.log("\n")
            self.log.header1("Save region file")
            self.log("\n")
            
            # Open file   
            f = open(self.ds9file,"w")
            
            # Write coordinate system
            f.write("fk5\n")
             
            # Loop over models
            for model in self.models:
                
                # Continue only if point source or extended source model
                if model.type() == "PointSource" or model.type() == "ExtendedSource":   
                    
                    line = self.model2ds9string(model)
                    if len(line):
                        f.write(line+"\n") 
                
                # Logging for diffuse components    
                elif model.type() == "DiffuseSource":
                    self.log("Skipping diffuse model "+model.name()+"\n")
                    
                # Logging for background components 
                else:
                    if self.logExplicit():
                        self.log("Skipping background model "+model.name()+"\n")   
                                      
            # Close file
            f.close()
        
        # Return
        return 

    def model2ds9string(self, model):
        """
        Converts model instance to ds9 region string 
        """
        
        # Initialise return string
        ds9string = ""
        
        # Determine color
        color = self.fixed_color
        for par in model:
            if par.is_free():
                color = self.free_color
                break
        
        # Initialise flag ot continue
        is_valid = True  
         
        # Retrieve model position  
        try:
            modelpos = model.spatial().dir()
        except AttributeError:            
            self.log("Skipping source "+model.name()+": Could not determine sky position\n")
            is_valid = False
        
        # Continue only if position was found
        if is_valid: 
            
            # Retrieve model type and name   
            modeltype = model.type()
            modelname = model.name()
            
            # Handle point source
            if modeltype == "PointSource":
                
                ds9string += "point("
                ds9string += str(modelpos.ra_deg()) 
                ds9string += ","
                ds9string += str(modelpos.dec_deg())
                ds9string += ") # "
                ds9string += "point="+self.pnt_type+" "
                ds9string += str(self.pnt_mark_size)+" "
            
            # Hanlde extended sources    
            elif modeltype == "ExtendedSource":
                
                # Retrieve spatial model
                spatial = model.spatial()
                classname = spatial.classname()
                
                # Handle radial sources
                if "Radial" in classname:  
                    
                    # Retrieve short name of model class (e.g. "Gauss" or "Disk)
                    shorttype = classname.split("Radial")[-1]
                    
                    # Handle Disk and Shell model
                    if classname == "GModelSpatialRadialDisk" or classname == "GModelSpatialRadialShell":
                        size = spatial.radius()
                    
                    # Handle Gauss Model
                    elif classname == "GModelSpatialRadialGauss":
                        size = spatial.sigma()
                        
                    # Skip if source is unknown
                    else:
                        self.log("Skipping source "+model.name()+": Spatial model \""+classname+"\" unkown\n")
                        is_valid = False
                
                    # Append radial ds9 string
                    if is_valid:
                        ds9string += "circle("
                        ds9string += str(modelpos.ra_deg()) 
                        ds9string += ","
                        ds9string += str(modelpos.dec_deg())
                        ds9string += ","
                        ds9string += str(size*3600.)+"\") #"          
                
                # Handle elliptical sources 
                elif "Elliptical" in classname:
                    
                    # Retrieve short name and source size
                    shorttype = classname.split("Elliptical")[-1]
                    size1 = spatial.semimajor()
                    size2 = spatial.semiminor()
                    angle = spatial.posangle()
    
                    # Append elliptical ds9 string
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
                    self.log("Skipping source "+model.name()+": Spatial model neither \"Radial\" nor \"Elliptical\"\n")   
                
                # Add short name to modelname
                if self.show_ext_type:
                    modelname +=" ("+shorttype+")"
    
            # Add user defined attributes
            ds9string += "color="+color+" "
            ds9string += "width="+str(self.width)+" "
            ds9string += "font=\""+self.font+" "+str(self.fontsize)+" "+self.fonttype+" "+self.fonttype2+"\" "
            if self.show_labels:
                ds9string += "text={"+modelname+"}"
        
        # Return string
        return ds9string 

    def execute(self):
        """
        Execute the script.
        """
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
    """
    Generates model information output.
    """
    # Create instance of application
    app = csmodelinfo(sys.argv)
    
    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
