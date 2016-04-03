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
    Shows the content of an observation container.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name    = "csobsinfo"
        self._version = "1.1.0"

        # Initialise class members
        self._obj_dir        = None
        self._compute_offset = False
        self._offsets        = []
        self._zeniths        = []
        self._azimuths       = []
        self._pnt_ra         = []
        self._pnt_dec        = []
        self._ebounds        = gammalib.GEbounds()
        self._gti            = gammalib.GGti()

        # Initialise observation container from constructor arguments.
        self._obs, argv = self._set_input_obs(argv)

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
        # Get parameters   
        if self._obs.size() == 0:  
            self._require_inobs("csobsinfo::get_parameters")
            self._obs = self._get_observations(False)
        
        # Initialise object position
        self._obj_dir = gammalib.GSkyDir()

        # Get (optional) offset parameters
        self._compute_offset = self["offset"].boolean()
        if self._compute_offset:
            ra  = self["ra"].real()
            dec = self["dec"].real() 
            self._obj_dir.radec_deg(ra,dec)

        # Read ahead DS9 filename
        if self._read_ahead():
            self["ds9file"].filename()

        # Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")

        # Return
        return
        

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
        
        # Initialise arrays to store certain values for reuse
        # Todo, think about using a python dictionary
        self._offsets  = []
        self._zeniths  = []
        self._azimuths = []
        self._pnt_ra   = []
        self._pnt_dec  = []
        self._ebounds  = gammalib.GEbounds()
        self._gti      = gammalib.GGti()
        obs_names      = []
        
        # Initialise output to be filled
        ontime         = 0.0
        livetime       = 0.0
        n_events       = 0
        n_eventbins    = 0
        n_obs_binned   = 0
        n_obs_unbinned = 0
        
        # Write header
        if self._logTerse():
            self._log("\n")
            if self._obs.size() > 1:
                self._log.header1("Observations")
            else:
                self._log.header1("Observation")
        
        # Loop over observations
        for obs in self._obs:

            # Skip non-CTA observations
            if not obs.classname() == "GCTAObservation":
                self._log("Skipping "+obs.instrument()+" observation\n")
                continue 
            
            # Use observed object as observation name if name is not given
            obs_name = obs.name()
            if obs_name == "":
                obs_name = obs.object()
            
            # Logging
            if self._logTerse():
                self._log.header2(obs_name)
                
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
            self._pnt_ra.append(pnt_dir.ra_deg())
            self._pnt_dec.append(pnt_dir.dec_deg())       
            
            # If avaliable append energy boundaries
            if obs_bounds.size() > 0 : 
                self._ebounds.append(obs_bounds.emin(),obs_bounds.emax())
            
            # Append time interval
            self._gti.append(obs_gti.tstart(), obs_gti.tstop())
            
            # Increment global livetime and ontime
            ontime   += obs.ontime()
            livetime += obs.livetime()
            
            # Bookkeeping
            if obs.eventtype() == "CountsCube":
                n_eventbins  += obs.events().size()
                n_obs_binned += 1
                if self._logTerse():
                    self._log.parformat("Binned")
                    self._log("yes\n")
                    self._log.parformat("Number of bins")
                    self._log(str(obs.events().size()))
                    self._log("\n")
            else:
                n_events       +=  obs.events().size()
                n_obs_unbinned += 1
                if self._logTerse():
                    self._log.parformat("Binned")
                    self._log("no\n")
                    self._log.parformat("Number of events")
                    self._log(str(obs.events().size()))
                    self._log("\n")

            # Retrieve zenith and azimuth and store for later use
            zenith  = obs.pointing().zenith()
            azimuth = obs.pointing().azimuth()
            self._zeniths.append(zenith)
            self._azimuths.append(azimuth)

            # Optionally compute offset with respect to target direction
            if self._compute_offset:
                offset = pnt_dir.dist_deg(self._obj_dir)
                self._offsets.append(offset)
            else:
                self._offsets.append(-1.0)
                    
            # Optionally log details
            if self._logTerse():
                
                # Log the observation energy range (if available)
                self._log.parformat("Energy range")
                if obs_bounds.size() == 0:
                    self._log("undefined")
                else:    
                    self._log(str(obs_bounds.emin()))
                    self._log(" - ")
                    self._log(str(obs_bounds.emax()))
                self._log("\n")
                
                # Log observation time interval
                self._log.parformat("Time range (MJD)")
                if obs_gti.size() == 0:
                    self._log("undefined")
                else:    
                    self._log(str(obs_gti.tstart().mjd()))
                    self._log(" - ")
                    self._log(str(obs_gti.tstop().mjd()))
                self._log("\n")
                
                # Log ontime
                self._log.parformat("Ontime")
                self._log(str(obs.ontime()))
                self._log(" s\n")
                
                # Log livetime
                self._log.parformat("Livetime")
                self._log(str(obs.livetime()))
                self._log(" s\n")
                
                # Log dead time fraction
                self._log.parformat("Deadtime fraction (%)")
                self._log("%.3f" % (deadfrac))
                self._log("\n")
                
                # Log pointing direction
                self._log.parformat("Pointing")
                self._log(str(pnt_dir))
                self._log("\n")
                
                # Optionally log offset with respect to target direction
                if self._compute_offset:
                    self._log.parformat("Offset from target")
                    self._log("%.2f" % (offset))
                    self._log(" deg\n")
                    
                # Log Zenith and Azimuth if required
                self._log.parformat("Zenith angle")
                self._log("%.2f" % (zenith))
                self._log(" deg\n")
                self._log.parformat("Azimuth angle")
                self._log("%.2f" % (azimuth))
                self._log(" deg\n")
                    
        # Log summary
        if self._logTerse():

            # Write header
            self._log("\n")
            self._log.header1("Summary")
        
            # Log general summary
            self._log.header3("Observations")
            self._log.parformat("Unbinned observations")
            self._log(str(n_obs_unbinned))
            self._log("\n")
            self._log.parformat("Binned observations")
            self._log(str(n_obs_binned))
            self._log("\n")
            self._log.header3("Events")
            self._log.parformat("Number of events")
            self._log(str(n_events))
            self._log("\n")
            self._log.parformat("Number of bins")
            self._log(str(n_eventbins))
            self._log("\n")
        
            # Log pointing summary
            self._log.header3("Pointings")
        
            # Log mean offset if possible
            if self._compute_offset:
                self._log.parformat("Mean offset angle")
                self._log("%.2f" % (sum(self._offsets)/len(self._offsets)))
                self._log(" deg\n")
        
            # Log mean azimuth and zenith angle
            self._log.parformat("Mean zenith angle")
            self._log("%.2f" % (sum(self._zeniths)/len(self._zeniths)))
            self._log(" deg\n")
            self._log.parformat("Mean azimuth angle")
            self._log("%.2f" % (sum(self._azimuths)/len(self._azimuths)))
            self._log(" deg\n")
        
            # Optionally log name of observations
            if self._logExplicit():
                obs_set = set(obs_names)
                for name in obs_set:
                    self._log.parformat("\""+name+"\"")
                    self._log(str(obs_names.count(name)))
                    self._log("\n")
                self._log("\n")
            
            # Log energy range
            self._log.header3("Energy range")
            self._log.parformat("Minimum energy")
            if self._ebounds.size() == 0:
                self._log("undefined")
            else:
                self._log(str(self._ebounds.emin()))
            self._log("\n")
            self._log.parformat("Maximum energy")
            if self._ebounds.size() == 0:
                self._log("undefined")
            else:
                self._log(str(self._ebounds.emax()))
            self._log("\n")
        
            # Log time range
            self._log.header3("Time range")
            self._log.parformat("Start (MJD)")
            self._log(str(self._gti.tstart().mjd()))
            self._log("\n")
            self._log.parformat("Stop (MJD)")
            self._log(str(self._gti.tstop().mjd()))
            self._log("\n")  
        
            # Log ontime and livetime in different units      
            self._log.parformat("Total ontime")
            self._log("%.2f s = %.2f min = %.2f h" % (ontime, ontime/60., ontime/3600.))
            self._log("\n")
            self._log.parformat("Total livetime")
            self._log("%.2f s = %.2f min = %.2f h" % (livetime, livetime/60.0, livetime/3600.))
            self._log("\n")        
                
        # Return
        return

    def save(self):
        """ 
        Save pointings into DS9 region file.

        This method saves all pointing directions that are found in the
        observation container into a DS9 region file. If "NONE" is
        specified for the "ds9file" parameter the method does nothing.
        """
        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Save pointings in DS9 file")

        # Get output filename in case it was not read ahead
        ds9file = self["ds9file"].filename()

        # Check if DS9 file is valid
        if ds9file.url() != "NONE":      

            # Log filename
            if self._logTerse():
                self._log(gammalib.parformat("DS9 filename"))
                self._log(ds9file.url())
                self._log("\n")
            
            # Open file   
            f = open(ds9file.url(),"w")
            
            # Write coordinate system
            f.write("fk5\n")
            
            # Loop over pointings
            for i in range(len(self._pnt_ra)):
                
                # Create string
                line = "point("
                line+=str(self._pnt_ra[i])+","+str(self._pnt_dec[i])+")"
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
        return self._zeniths
    
    def azimuths(self):
        """
        Return azimuth angles
        """
        return self._azimuths
    
    def offsets(self):
        """
        Return offset angles
        """
        return self._offsets

    def ebounds(self):
        """
        Return energy boundaries
        """
        return self._ebounds
    
    def gti(self):
        """
        Return good time intervals
        """
        return self._gti
            
    def execute(self):
        """
        Execute the script.
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
    app = csobsinfo(sys.argv)
    
    # Execute application
    app.execute()
