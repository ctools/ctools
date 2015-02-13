#! /usr/bin/env python
# ==========================================================================
# spectral points generation script.
#
# Copyright (C) 2014-2015 Michael Mayer
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

try:
    import matplotlib.pyplot as plt
except:
    sys.exit("Need matplotlib to launch this script")

# Initialise conversion constant
MeV2erg = 1.60217e-6

# ============== #
# cstsdist class #
# ============== #
class csplotspec(ctools.cscript):
    """
    This class plots a copmplete spectrum. It derives from
    the ctools.cscript class which provides support for parameter files,
    command line arguments, and logging. In that way the Python script
    behaves just as a regular ctool. 
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        
        # Set name
        self.name    = "csplotspec"
        self.version = "1.0.0"
        self.outfile = ""
        
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
        #  Write separator into logger
        if self.logTerse():
            self.log("\n")
        
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
            pars.append(gammalib.GApplicationPar("inmodel","f","a","results.xml","","","Source model"))
            pars.append(gammalib.GApplicationPar("srcname","s","a","Crab","","","Source name"))
            pars.append(gammalib.GApplicationPar("inspec","f","a","spectrum.fits","","","Input spectral points"))
            pars.append(gammalib.GApplicationPar("inbfly","f","a","butterfly.txt","","","Input butterfly"))
            pars.append(gammalib.GApplicationPar("outfile","f","a","spectrum.png","","","Output file name"))
            pars.append(gammalib.GApplicationPar("tslim","r","h","4.0","","","Minimum TS value for data point"))
            pars.append(gammalib.GApplicationPar("xmin","r","h","0.1","","","Minimum energy for plotting"))
            pars.append(gammalib.GApplicationPar("xmax","r","h","100","","","Maximum energy for plotting"))
            pars.append(gammalib.GApplicationPar("show","b","h","yes","yes/no","","Show plot directly"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csplotspec.log","","","Log filename"))
            pars.save(parfile)
        
        # Return
        return
        
    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get parameters     
        self.models = gammalib.GModels(self["inmodel"].filename())
        self.m_source = self.models[self["srcname"].string()]
        self.m_bfly = self["inbfly"].filename()
        self.m_inspec = self["inspec"].filename()
        self.m_outfile = self["outfile"].filename()            
             
        # Get other parameters     
        self.m_tslim = self["tslim"].real()
        self.m_xmin = self["xmin"].real()
        self.m_xmax = self["xmax"].real()
        self.m_show = self["show"].boolean()
             
        # Set some fixed parameters
        self.m_log     = False # Logging in client tools
        self.m_chatter = self["chatter"].integer()
        self.m_clobber = self["clobber"].boolean()
        self.m_debug   = self["debug"].boolean()
         
        # Return
        return
    
        
    def get_points(self):
        
        # Get FITS file
        fits = gammalib.GFits(self.m_inspec)
        table = fits["SPECTRUM"]
        
        # Get FITS columns from file
        c_energy = table["energy"]
        c_flux = table["flux"]
        c_flux_err = table["flux_err"]
        c_TS = table["TS"]
        c_ulim = table["ulimit"]        
        
        # Initialise arrays
        energy = []
        flux = []
        ed_flux = []
        eu_flux = []
        
        # Upper limit arrays
        lim_e = []
        lim_flux = []
        lim_edflux = []
        lim_euflux = []
        
        # Loop over points
        for i in range(c_energy.length()):
                
            # Check for significant point
            if c_TS[i] < self.m_tslim:

                # Use upper limit if existent
                if c_ulim[i] > 0.0:
                    
                    # Append value to array
                    lim_e.append(c_energy[i])
                    lim_flux.append(c_ulim[i])
                    lim_edflux.append(0.9*c_ulim[i])
                    lim_euflux.append(0.0)
                else:
                    # If there is no upper limit available, this point is skipped
                    continue
                    
            else:
                
                # Append significant data point to array
                energy.append(c_energy[i])
                flux.append(c_flux[i])
                ed_flux.append(c_flux_err[i])
                eu_flux.append(c_flux_err[i])
 
        # Return plot values
        return energy,flux,ed_flux,eu_flux, lim_e, lim_flux, lim_edflux, lim_euflux

    def get_line(self, xmin = 0.1, xmax = 100.):
        
        # Initialise energy binning
        ebounds = gammalib.GEbounds(100, gammalib.GEnergy(xmin,"TeV"), gammalib.GEnergy(xmax,"TeV"), True)
        
        # Initialise return array 
        x = []
        y = []
        
        # Initialise dummy time
        t = gammalib.GTime()
        
        # Loop over energy bins
        for i in range(ebounds.size()):
            
            # Append values to array
            energy = ebounds.elogmean(i)
            x.append(energy.TeV())
            y.append(energy.MeV() * energy.MeV() * self.m_source.spectral().eval(energy, t) * MeV2erg)
        
        # Return plot values
        return x,y
     
    def get_butterfly(self):
        
        # Get ascii file
        csv = gammalib.GCsv(self.m_bfly)
        
        # Initialise return arrays
        x = []
        y = []
        
        # Loop over rows
        for i in range(csv.nrows()):
            
            #store energy
            energy = csv.real(i,0) # MeV
            
            # store upper limit flux
            flux = csv.real(i,1) + csv.real(i,2) / 2.0 #ph/cm2/s/MeV 
            
            # append to array           
            x.append(energy * 1e-6) # TeV
            y.append(flux * energy * energy * MeV2erg) # erg/cm2/s
            
        # Loop again over file from the other side   
        for i in range(csv.nrows()):
            
            # Create index counting backwards
            k = csv.nrows()-1-i
            
            # Get values
            energy = csv.real(k,0) # MeV
            flux = csv.real(k,1)-csv.real(k,2) / 2.0 #ph/cm2/s/MeV    
            
            # Append values        
            x.append(energy * 1e-6) # TeV
            y.append(flux * energy * energy * MeV2erg) # erg/cm2/s
        
        # Return plot values
        return x,y
        
               
        
    def execute(self):
        """
        Execute the script.
        """
        # Run the script
        self.run()

        # Save image file
        if not (self.m_outfile == "" or self.m_outfile == "NONE"):
            plt.savefig(self.outfile)
        
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

        plt.figure("Spectrum")
        plt.loglog()
        plt.xlabel("Energy [TeV]")
        plt.ylabel(r'$\nu$F$\nu$ [erg cm$^{-2}$ s$^{-1}$]')

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Get Butterfly")

        if not (self.m_bfly == "NONE" or self.m_bfly == ""): 
            xbf, ybf = self.get_butterfly()   
            plt.fill(xbf,ybf,color='green',alpha=0.5)       

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Get Spectrum Line")
        
        xline,yline = self.get_line(self.m_xmin, self.m_xmax)    
        plt.plot(xline,yline, color='black',lw=2.0)   

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Get Spectrum Points")
        
        if not (self.m_inspec == "NONE" or self.m_inspec == ""):     
            e,f,edf,euf,le,lf,ledf,leuf = self.get_points()   
            plt.errorbar(x=e,y=f,yerr=[edf,euf],fmt='bo',lw=2.0)
            if len(le):
                plt.errorbar(x=le,y=lf,yerr=[ledf,leuf],fmt='bo',lw=2.0, capsize=2.0,lowlims=True)

        if self.m_show:
            plt.show()

            
        # Return
        return
    

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Generates spectral points.
    """
    # Create instance of application
    app = csplotspec(sys.argv)
    
    # Open logfile
    app.logFileOpen()
    
    # Execute application
    app.execute()
    