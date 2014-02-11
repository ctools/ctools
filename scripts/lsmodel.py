#! /usr/bin/env python
# ==========================================================================
# This script generates a model map from a counts map and a source model.
#
# Copyright (C) 2014 Juergen Knoedlseder
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
import sys


# ============= #
# lsmodel class #
# ============= #
class lsmodel(gammalib.GApplication):
    """
    This class implements the model map generation script. It derives from
    the GammaLib::GApplication class which provides support for parameter
    files, command line arguments, and logging. In that way the Python
    script behaves just as a regular ctool.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "lsmodel"
        self.version = "0.1.0"
        
        # Make sure that parfile exists
        file = self.parfile()

        # Initialise application
        if len(argv) == 0:
            gammalib.GApplication.__init__(self, self.name, self.version)
        elif len(argv) ==1:
            gammalib.GApplication.__init__(self, self.name, self.version, *argv)
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
        #  Write separator into logger
        if self.logTerse():
            self.log("\n")
        
        # Return
        return

    def parfile(self):
        """
        Check if parfile exists. If parfile does not exist then create a
        default parfile. This kluge avoids shipping the script with a parfile.
        """
        # Set parfile name
        parfile = self.name+".par"
        
        try:
            pars = gammalib.GApplicationPars(parfile)
        except:
            # Signal if parfile was not found
            sys.stdout.write("Parfile "+parfile+" not found. Create default parfile.\n")
            
            # Create default parfile
            pars = gammalib.GApplicationPars()
            pars.append(gammalib.GApplicationPar("cntmap","f","a","cntmap.fits","","","Counts map"))
            pars.append(gammalib.GApplicationPar("expmap","f","a","binned_expmap.fits","","","Exposure map"))
            pars.append(gammalib.GApplicationPar("ltcube","f","a","ltcube.fits","","","Livetime cube"))
            pars.append(gammalib.GApplicationPar("srcmdl","f","a","crab_model.xml","","","Source model"))
            pars.append(gammalib.GApplicationPar("caldb","s","h","","","","Calibration database"))
            pars.append(gammalib.GApplicationPar("irf","s","a","P7REP_SOURCE_V15","","","Instrument response function"))
            pars.append(gammalib.GApplicationPar("outfile","f","a","model_cube.fits","","","Output file name"))
            pars.append_standard()
            pars.save(parfile)
        
        # Return
        return
        
    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get parameters
        self.m_cntmap  = self["cntmap"].filename()
        self.m_expmap  = self["expmap"].filename()
        self.m_ltcube  = self["ltcube"].filename()
        self.m_srcmdl  = self["srcmdl"].filename()
        self.m_caldb   = self["caldb"].string()
        self.m_irf     = self["irf"].string()
        self.m_outfile = self["outfile"].filename()

        # Set some fixed parameters
        self.m_log   = False # Logging in client tools
        self.m_debug = False # Debugging in client tools
        
        # Return
        return
    
    def execute(self):
        """
        Execute the script.
        """
        # Run the script
        self.run()
        
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
        
        #  Write input parameters into logger
        if self.logTerse():
            self.log_parameters()
            self.log("\n")
        
        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Generate model cube")

        # Set LAT observation
        obs = gammalib.GLATObservation()
        obs.load_binned(self.m_cntmap, self.m_expmap, self.m_ltcube)
        obs.response(self.m_irf, self.m_caldb)

        # Set model (use first model for the moment)
        models = gammalib.GModels(self.m_srcmdl)
        model  = models[0]

        # Get deep copy of GLATEventCube from observation
        cube = gammalib.GLATEventCube(obs.events())

        # Fill cube with model values
        for i in range(cube.size()):
            event = cube[i]
            value = model.eval(event, obs)
            event.counts(value)
            cube[i] = event

        # Log cube
        if self.logTerse():
            self.log(str(cube))
            self.log("\n")

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Save model cube")

        # Save model cube
        cube.save(self.m_outfile, self.clobber())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Generates a model map from a counts map and a source model.
    """
    # Create instance of application
    app = lsmodel(sys.argv)
    
    # Open logfile
    app.logFileOpen()
    
    # Execute application
    app.execute()
    