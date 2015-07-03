#! /usr/bin/env python
# ==========================================================================
# This script dumps all available calibrations into the console.
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
import gammalib
import ctools
import sys
import glob
import os


# ============ #
# cspull class #
# ============ #
class cscaldb(ctools.cscript):
    """
    This class dumps the content of the ctools calibration database.
    It derives from the ctools.cscript class which provides support for
    parameter files, command line arguments, and logging. In that way
    the Python script behaves just as a regular ctool.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "cscaldb"
        self.version = "1.0.0"
        
        # Make sure that parfile exists
        file = self.parfile()

        # Initialise application
        if len(argv) == 0:
            ctools.cscript.__init__(self, self.name, self.version)
        elif len(argv) ==1:
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
        default parfile. This kluge avoids shipping the cscript with a parfile.
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
            pars.append_standard()
            pars.save(parfile)
        
        # Return
        return
        
    def get_parameters(self):
        """
        Get parameters from parfile.
        """
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

        # Get the calibration database path
        caldb = gammalib.GCaldb()

        # Gather missions
        missions = []
        paths    = glob.glob(caldb.rootdir()+"/data/*")
        paths.sort()
        for path in paths:
            head, tail = os.path.split(path)
            missions.append(tail)

        # Loop over missions
        for mission in missions:

            # Skip all non-CTA instruments
            if mission != "cta":
                continue

            # Write mission into logger
            if self.logTerse():
                self.log("\n")
                self.log.header1("Mission: "+mission)

            # Gather instruments
            instruments = []
            paths       = glob.glob(caldb.rootdir()+"/data/"+mission+"/*")
            paths.sort()
            for path in paths:
                head, tail = os.path.split(path)
                instruments.append(tail)

            # Loop over instruments
            for instrument in instruments:
            
                # Write mission into logger
                if self.logTerse():
                    self.log.header3("Response functions in database \""+instrument+"\"")

                # Open calibration index file
                cifname = caldb.rootdir()+"/data/"+mission+"/"+instrument+"/caldb.indx"
                fits = gammalib.GFits(cifname)
                cif  = fits["CIF"]
                cals = cif["CAL_CBD"]
                
                # Extract response names
                names = []
                nrows = cals.length()
                ncols = cals.number()
                for row in range(nrows):
                    for col in range(ncols):
                        cal    = cals.string(row, col)
                        istart = cal.find("NAME(")
                        if istart != -1:
                            istop = cal.find(")")
                            name  = cal[5:istop]
                            if names.count(name) == 0:
                                names.append(name)
                names.sort()
                
                # Print response name
                if self.logTerse():
                    for name in names:
                        self.log(name+"\n")
                    self.log("\n")
        
        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Dumps calibration database into log file.
    """
    # Create instance of application
    app = cscaldb(sys.argv)
    
    # Open logfile
    app.logFileOpen()
    
    # Execute application
    app.execute()
    