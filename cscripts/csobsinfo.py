#!/usr/bin/env python
# ==========================================================================
# Dump information about observation into log file
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


# =============== #
# csobsinfo class #
# =============== #
class csobsinfo(ctools.cscript):
    """
    This class dumps information about an observation container into a
    logfile or on screen. This might be helpful for quick access to an
    observation container showing, e.g., total lifetime, energy range etc.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name    = "csobsinfo"
        self.version = "1.1.0"

        # Initialise some members
        if len(argv) > 0 and isinstance(argv[0],gammalib.GObservations):
            self.obs = argv[0]
            argv     = argv[1:]
        else:      
            self.obs = gammalib.GObservations()
            self.obs.clear()   

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
            pars.append(gammalib.GApplicationPar("inobs","f","a","obs.xml","","","Event list, counts cube, or observation definition file"))
            pars.append(gammalib.GApplicationPar("offset","b","h","no","","","Compute offset from target to pointing positions"))
            pars.append(gammalib.GApplicationPar("ra","r","a","0.0","","","Target Right Ascension"))
            pars.append(gammalib.GApplicationPar("dec","r","a","0.0","","","Target Declination"))
            pars.append(gammalib.GApplicationPar("ds9file","f","h","NONE","","","DS9 region file containing pointings"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csobsinfo.log","","","Log filename"))
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get parameters   
        if self.obs.size() == 0:  
            self.require_inobs("csobsinfo::get_parameters")
            self.obs = self.get_observations(False)
        
        # Initialise object position
        self.obj_dir = gammalib.GSkyDir()

        # Get (optional) offset parameters
        self.m_offset = self["offset"].boolean()
        if self.m_offset:
            ra  = self["ra"].real()
            dec = self["dec"].real() 
            self.obj_dir.radec_deg(ra,dec)

        # Get (optional) DS9 filename
        self.ds9file = self["ds9file"].filename()

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
        
        # Initialise arrays to store certain values for reuse
        # Todo, think about using a python dictionary
        self.m_offsets  = []
        self.m_zeniths  = []
        self.m_azimuths = []
        self.pnt_ra   = []
        self.pnt_dec  = []
        obs_names     = []
        
        # Initialise output to be filled
        self.m_ebounds   = gammalib.GEbounds()
        self.m_gti       = gammalib.GGti()
        ontime         = 0.0
        livetime       = 0.0
        n_events       = 0
        n_eventbins    = 0
        n_obs_binned   = 0
        n_obs_unbinned = 0
        
        # Logging
        self.log.header1("Looping over observations")
        
        # Loop over observations
        for obs in self.obs:
            

            # Skip non-CTA observations
            if not obs.classname() == "GCTAObservation":
                self.log("Skipping "+obs.instrument()+" observation\n")
                continue 
            
            # Use observed object as observation name
            # if name is not given
            obs_name = obs.name()
            if obs_name == "":
                obs_name = obs.object()
            
            # Logging
            if self.logTerse():
                self.log.header2(obs_name)
                
            # Retrieve observation name
            obs_names.append(obs_name)
            
            # Retrieve energy boundaries
            obs_bounds = obs.events().ebounds()
            
            # Retrieve time interval
            obs_gti = obs.events().gti()
            
            # Compute mean time and dead time fraction in percent
            deadfrac = (1.0-obs.deadc())*100.0
            
            # Retrieve pointing and store Ra,Dec
            pnt_dir = obs.pointing().dir()  
            self.pnt_ra.append(pnt_dir.ra_deg())
            self.pnt_dec.append(pnt_dir.dec_deg())       
            
            # If avaliable append energy boundaries
            if obs_bounds.size() > 0 : 
                self.m_ebounds.append(obs_bounds.emin(),obs_bounds.emax())
            
            # Append time interval
            self.m_gti.append(obs_gti.tstart(), obs_gti.tstop())
            
            # Increment global livetime and ontime
            ontime   += obs.ontime()
            livetime += obs.livetime()
            
            # Log the event type (binned or unbinned)
            if self.logTerse():
                self.log.parformat("binned")
            
            if obs.eventtype() == "CountsCube":
                n_eventbins  += obs.events().size()
                n_obs_binned += 1
                if self.logTerse():
                    self.log("yes")
                    self.log("\n")
                    self.log.parformat("Number of bins")
                    self.log(str(obs.events().size()))
                    self.log("\n")
            else:
                n_events       +=  obs.events().size()
                n_obs_unbinned += 1
                if self.logTerse():
                    self.log("no")
                    self.log("\n")
                    self.log.parformat("Number of events")
                    self.log(str(obs.events().size()))
                    self.log("\n")
                    
            # Log option
            if self.logTerse():  
                
                # Log the observation energy range (if available)
                self.log.parformat("Energy range")
                if obs_bounds.size() == 0:
                    self.log("undefined")
                else:    
                    self.log(str(obs_bounds.emin())+" - "+str(obs_bounds.emax()))
                self.log("\n")
                
                # Log observation time interval
                self.log.parformat("Time range (MJD)")
                if obs_gti.size() == 0:
                    self.log("undefined")
                else:    
                    self.log(str(obs_gti.tstart().mjd())+" - "+str(obs_gti.tstop().mjd()))
                self.log("\n")
                
                # Log ontime
                self.log.parformat("Ontime [s]")
                self.log(str(obs.ontime()))
                self.log("\n")
                
                # Log livetime
                self.log.parformat("Livetime [s]")
                self.log(str(obs.livetime()))
                self.log("\n")
                
                # Log dead time fraction
                self.log.parformat("Dead-time fraction (%)")
                self.log("%.3f" % (deadfrac))
                self.log("\n")
                
                # Log pointing direction
                self.log.parformat("Pointing")
                self.log(str(pnt_dir))
                self.log("\n")
                
                # Log offset if possible
                if self.m_offset:
                    offset = pnt_dir.dist_deg(self.obj_dir)
                    self.m_offsets.append(offset)
                    self.log.parformat("Offset")
                    self.log("%.2f" % (offset))
                    self.log("\n")
                else:
                    self.m_offsets.append(-1.0)
                    
                # Retrieve zenith and azimuth and store for later use
                zenith = obs.pointing().zenith()
                azimuth = obs.pointing().azimuth()
                self.m_zeniths.append(zenith)
                self.m_azimuths.append(azimuth)
                
                # Log Zenith and Azimuth if required
                if self.logExplicit():
                    self.log.parformat("Zenith angle")
                    self.log("%.2f" % (zenith))
                    self.log("\n")
                    self.log.parformat("Azimuth angle")
                    self.log("%.2f" % (azimuth))
                    self.log("\n\n")
                    
        # Log summary
        self.log("\n")
        self.log.header1("Summary")
        
        # Log general summary
        self.log.parformat("Unbinned observations")
        self.log(str(n_obs_unbinned))
        self.log("\n")
        self.log.parformat("Number of events")
        self.log(str(n_events))
        self.log("\n")
        self.log.parformat("Binned observations")
        self.log(str(n_obs_binned))
        self.log("\n")
        self.log.parformat("Number of bins")
        self.log(str(n_eventbins))
        self.log("\n\n")
        
        # Log pointing summary
        self.log.header3("Pointings")
        
        # Log mean offset if possible
        if self.m_offset:
            self.log.parformat("Mean offset angle")
            self.log("%.2f" % (sum(self.m_offsets)/len(self.m_offsets)))
            self.log("\n")
        
        # Log mean azimuth and zenith angle
        self.log.parformat("Mean zenith angle")
        self.log("%.2f" % (sum(self.m_zeniths)/len(self.m_zeniths)))
        self.log("\n")
        self.log.parformat("Mean azimuth angle")
        self.log("%.2f" % (sum(self.m_azimuths)/len(self.m_azimuths)))
        self.log("\n")
        
        # Log name of observations if requested
        if self.logExplicit():
            obs_set = set(obs_names)
            for name in obs_set:
                self.log.parformat("\""+name+"\"")
                self.log(str(obs_names.count(name)))
                self.log("\n")
        self.log("\n")
            
        # Log energy range
        self.log.header3("Energy range")
        self.log.parformat("Emin")
        if self.m_ebounds.size() == 0:
            self.log("undefined")
        else:
            self.log(str(self.m_ebounds.emin()))
        self.log("\n")
        self.log.parformat("Emax")
        if self.m_ebounds.size() == 0:
            self.log("undefined")
        else:
            self.log(str(self.m_ebounds.emax()))
        self.log("\n\n")
        
        # Log time range
        self.log.header3("Time range")
        self.log.parformat("Start [MJD]")
        self.log(str(self.m_gti.tstart().mjd()))
        self.log("\n")
        self.log.parformat("Stop [MJD]")
        self.log(str(self.m_gti.tstop().mjd()))
        self.log("\n")  
        
        # Log ontime and livetime in different units      
        self.log.parformat("Total ontime")
        self.log("%.2f s = %.2f min = %.2f h" % (ontime, ontime/60., ontime/3600.))
        self.log("\n")
        self.log.parformat("Total livetime")
        self.log("%.2f s = %.2f min = %.2f h" % (livetime, livetime/60.0, livetime/3600.))
        self.log("\n")        
                
        # Return
        return

    def save(self):
        """ 
        Save pointings to ds9 region file if required
        """
        
        # Check if ds9 file is valid
        if not self.ds9file == "NONE":      
            
            # Open file   
            f = open(self.ds9file.url(),"w")
            
            # Write coordinate system
            f.write("fk5\n")
            
            # Loop over pointings
            for i in range(len(self.pnt_ra)):
                
                # Create string
                line = "point("
                line+=str(self.pnt_ra[i])+","+str(self.pnt_dec[i])+")"
                line+=" # point=cross 20 width=3 \n"
                
                # Write to file
                f.write(line)
                
            # Close file
            f.close()
        
        # Return
        return 
            
    
    def zeniths(self):
        """
        Return zenith angles
        """
        return self.m_zeniths
    
    def azimuths(self):
        """
        Return azimuth angles
        """
        return self.m_azimuths
    
    def offsets(self):
        """
        Return offset angles
        """
        return self.m_offsets

    def ebounds(self):
        """
        Return energy boundaries
        """
        return self.m_ebounds
    
    def gti(self):
        """
        Return good time intervals
        """
        return self.m_gti
            
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
    Generates information about an observation container.
    """
    # Create instance of application
    app = csobsinfo(sys.argv)
    
    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
