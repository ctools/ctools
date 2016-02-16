#!/usr/bin/env python
# ==========================================================================
# Download IACT data from remote machine
#
# Copyright (C) 2016 Michael Mayer
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
import os
import json
import shutil

# ================ #
# csiactdata class #
# ================ #
class csiactcopy(ctools.cscript):
    """
    This script copies IACT data from one location to another.
    It can take a list of observation IDs to allow the download specific
    observations. Index files get merged and updated accordingly.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self.name     = "csiactcopy"
        self.version  = "1.1.0"
        self.datapath = os.getenv("VHEFITS","")

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
            # Signal if parfile was not found
            sys.stdout.write("Parfile "+parfile+" not found. Create default parfile.\n")

            # Create default parfile
            pars = gammalib.GApplicationPars()    
            pars.append(gammalib.GApplicationPar("remote_master","f","a","/Volumes/mountpoint/master.json","","","Location of remote master file"))      
            pars.append(gammalib.GApplicationPar("prodname","s","a","iact-fits","","","Name of FITS production to download"))   
            pars.append(gammalib.GApplicationPar("outpath","s","a","/path/to/fits","","","Destination path of FITS data"))        
            pars.append(gammalib.GApplicationPar("runlist","f","h","NONE","","","List of observation IDs"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csiactcopy.log","","","Log filename"))
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        
        # Get Parameters
        self.m_remote_master = self["remote_master"].filename()
        if not self.m_remote_master.exists():
            raise RuntimeError("*** ERROR: FITS data not available. No master index file found in \""+self.m_remote_master.url()+"\". Make sure remote file system is properly mounted")
        
        # Get parameters
        self.m_prodname = self["prodname"].string()
        self.m_outpath  = gammalib.expand_env(self["outpath"].string())
        self.m_runlist  = self["runlist"].filename()            
        
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
        # Switch screen logging on
        self.log.cout(True)

        # Get parameters
        self.get_parameters()

        # Write input parameters into logger
        if self.logTerse():
            self.log_parameters()
            self.log("\n")
        
        # Make destination directory if not available
        if not os.path.isdir(self.m_outpath):            
            os.makedirs(self.m_outpath)   
        
        # Log output
        if self.logTerse():
            self.log("\n")
            self.log.header2("Runlist information")
        
        # Initialise runlist array   
        self.runlist = []
        
        # Check if runlist file is available
        if self.m_runlist.exists():
            runfile = open(self.m_runlist)
            for line in runfile.readlines():
                if len(line) == 0:
                    continue
                if line[0] == '#':
                    continue
                if len(line.split()) > 0:
                    self.runlist.append(int(line.split()[0]))
            runfile.close()   
            
            # Logging
            if self.logTerse():
                self.log.parformat("Number of observations")
                self.log(str(len(self.runlist)))
                self.log("\n")
        
        # Copy all data if runlist file is "NONE"
        elif self.m_runlist == "NONE":
            if self.logTerse():
                self.log("Copy all available data")
                self.log("\n")
        
        # Raise exception if file not valid
        else:
            raise RuntimeError("*** ERROR: Runlist file \""+self.m_runlist+"\" not available")
        
        # Check for availability of remote master file
        if not self.m_remote_master.exists():
            raise RuntimeError("*** ERROR: Remote master file \""+self.m_remote_master.url()+"\" does not exist")
        
        # Retrieve json data from remote master
        json_data = open(self.m_remote_master.url()).read()
        data      = json.loads(json_data) 
        if not "datasets" in data:
            raise RuntimeError("*** ERROR: Key \"datasets\" not available in remote master index file \""+self.m_remote_master.url()+"\".")
        
        # Get array of configurations
        configs   = data["datasets"]
        
        # Get remote paths
        self.remote_base = os.path.dirname(self.m_remote_master.url())
                
        # Initialise flag if prod has been found
        has_prod = False
        
        # Initialise file names to be copied
        files = set()
        
        # Logging
        if self.logTerse():
            self.log("\n")
            self.log.header2("Loop over remote configs")
        
        # Loop over configs
        for config in configs: 
            
            # Logging          
            if self.logVerbose():
                self.log.header2(str(config["name"]))
            
            # Check if prodname was found
            if config["name"] == self.m_prodname:
                
                # Inidicate that prodname has been found
                has_prod = True           
                
                # Build path of index files
                remote_hdu = os.path.join(self.remote_base, config["hduindx"])
                remote_obs = os.path.join(self.remote_base, config["obsindx"])
                
                # Log information
                if self.logTerse():
                    self.log.header3("Remote config \""+self.m_prodname+"\"")
                    self.log.parformat("HDU index")
                    self.log(str(remote_hdu))
                    self.log("\n")
                    self.log.parformat("Observation index")
                    self.log(str(remote_obs))
                    self.log("\n")
                    
                # Open remote HDU index file
                fits = gammalib.GFits(str(remote_hdu))
                table = fits["HDU_INDEX"]
                
                # Initialise flag if SIZE column is present
                has_size = table.contains("SIZE")
                
                if has_size:
                    # Initialise file size
                    cp_size = 0
                
                # Initialise remote observation IDs
                remote_ids = set()
                for row in range(table.nrows()):
                    remote_ids.add(table["OBS_ID"][row])
                
                # Log runs that are not available remotely
                for run in self.runlist:
                    if not run in remote_ids:  
                        if self.logNormal():
                            self.log("Skip observation "+str(run)+": ID not available remotely")
                            self.log("\n")
                          
                
                # Loop over remote HDU index file
                for row in range(table.nrows()):
                    
                    # Get observation ID
                    obs_id = table["OBS_ID"][row]   
                    file_dir = table["FILE_DIR"][row]              
                    file_name = table["FILE_NAME"][row] 
                    
                    # Skip if filename is empty
                    if file_name == "":
                        continue
                                        
                    # Check if we need to consider an input runlist
                    if len(self.runlist):
                        
                        # Check if obs id is in runlist
                        if obs_id in self.runlist:                                    
                            
                            # Get filename
                            fname = os.path.join(os.path.dirname(remote_hdu), file_dir, file_name)
                            
                            # Add file
                            oldlen = len(files)
                            files.add(fname)
                            newlen = len(files)
                            
                            # Add file size
                            if has_size and newlen > oldlen:
                                cp_size += table["SIZE"][row]
                        
                    # Otherwise add every file       
                    else:
                        # Get filename
                        fname = os.path.join(os.path.dirname(remote_hdu), file_dir, file_name)
                        
                        # Add file
                        oldlen = len(files)
                        files.add(fname)
                        newlen = len(files)
                        
                        # Add file size
                        if has_size and newlen > oldlen:
                            cp_size += table["SIZE"][row]
                            
                # Log file information
                if self.logTerse():
                    self.log("\n")
                    self.log.header2("File information")
                    self.log.parformat("Number of files")
                    self.log(str(len(files)))
                    self.log("\n")
                    if has_size:
                        self.log.parformat("Size")
                        self.log("%.2f"%(float(cp_size)*1e-6)+" MB")
                        self.log("\n")
                if self.logVerbose():
                    self.log("\n")
                    self.log.header3("File names")
                    for filename in files:
                        self.log(str(filename)+"\n")
                
                # Close HDU index file         
                fits.close()
                
            # If prodname is not found just log that we skip the config
            else:
                if self.logExplicit():
                    self.log("\n")
                    self.log.header3("Skipping config \""+str(config["name"])+"\"")
                    self.log("\n")
                    
        # Raise Exception if prodname was not found                                    
        if not has_prod: 
            msg = "*** ERROR: FITS production \""+self.m_prodname+"\" not available. "
            msg += "Available productions are: \n"
            for config in configs:
                msg += " - "+config["name"] + "\n"
            raise RuntimeError(msg)
                
        # Logging
        if self.logNormal():
            self.log("\n")
            self.log.header1("Copying files")
        
        # intialise counter
        k = 0
    
        # initialise values for logging
        last_fraction = 0.0
        fraction_increment = 20.0
        if self.logNormal():
            fraction_increment = 10.0
        if self.logTerse():
            fraction_increment = 5.0
        if self.logExplicit():
            fraction_increment = 2.0
        
        # initialise logging properties
        n_copied = 0
        total_size = 0.0
        
        # Loop over files and copy
        for filename in files:
            
            # Log progress
            fraction = float(k) / float(len(files)) * 100.0
            while fraction > last_fraction:
                if self.logNormal() and not self.logVerbose():
                    self.log.parformat("Status")
                    self.log(str(int(last_fraction))+"%")
                    self.log("\n")
                last_fraction += fraction_increment
            
            # Copy file
            filesize = self.copy(filename, self.clobber())
            
            if filesize > 0.0:
                total_size += filesize
                n_copied   += 1
            
            # Increment counter
            k += 1
        
        # Logging
        if self.logNormal() and not self.logVerbose():
            self.log.parformat("Status")
            self.log("Finished")
            self.log("\n")
        if self.logTerse():
            self.log("\n")
            self.log.header1("Updating index files")
        
        # Build local index file names
        local_hdu = os.path.join(self.m_outpath, os.path.relpath(remote_hdu, self.remote_base))
        local_obs = os.path.join(self.m_outpath, os.path.relpath(remote_obs, self.remote_base))
        
        # If we have a runlist then merge index file
        if len(self.runlist):
            
            # Logging            
            if self.logTerse():
                self.log("\n")
                self.log.header3("HDU index")
                
            # Merge remote index files with local files
            self.merge(local_hdu, remote_hdu, "HDU_INDEX", self.clobber())
            
            if self.logTerse():
                self.log("\n")
                self.log.header3("OBS index")    
                
            # Merge remote index files with local files
            self.merge(local_obs, remote_obs, "OBS_INDEX", self.clobber())
            
        else: 
            # If all files were copied, just copy index files too
            self.copy(remote_hdu, self.clobber())
            self.copy(remote_obs, self.clobber())
  
        if self.logTerse():
            self.log("\n")
            self.log.header3("Master index file")
        
        # Adding prodname to local master
        localmaster = os.path.join(self.m_outpath, "master.json")
        
        # If no local master is found, copy master over first
        if not os.path.isfile(localmaster):
            self.copy(self.m_remote_master.url(), self.clobber())
            
        # Load local master
        json_data = open(localmaster).read()
        data      = json.loads(json_data) 
        configs   = data["datasets"]
        
        # Initialise flag indicating if we already have prodname in master file
        has_config = False
        
        # Initialise new configs array
        newconfigs = []
        
        # Loop over configs in master index file
        for config in configs:
            
            # Get hdu and obs index files
            hdu = os.path.join(self.m_outpath, config["hduindx"])
            obs = os.path.join(self.m_outpath, config["obsindx"])
            
            # Check if index files are available
            if not (gammalib.GFilename(str(hdu)).is_fits() and
                    gammalib.GFilename(str(obs)).is_fits()):
                if self.logTerse():
                    self.log("Removing \""+str(config["name"])+"\" (not available)")
                    self.log("\n")
            else:
                # Append config if valid
                newconfigs.append(config)
                if self.logTerse():
                    self.log("Keeping \""+str(config["name"])+"\"")
                    self.log("\n")
            
            # Signals that downloaded config is available          
            if config["name"] == self.m_prodname:
                has_config = True
        
        # Create new entry if config was not available    
        if not has_config:
            newdict = dict.fromkeys(["name","hduindx","obsindx"])
            newdict["name"] = self.m_prodname
            newdict["hduindx"] = os.path.relpath(local_hdu, self.m_outpath)
            newdict["obsindx"] = os.path.relpath(local_obs, self.m_outpath)
            newconfigs.append(newdict)
            if self.logTerse():
                self.log("Adding \""+str(newdict["name"])+"\"")
                self.log("\n")
        
        # Write new json master file
        f = open(localmaster, "w")
        data["datasets"] = newconfigs
        json.dump(data, f, indent=2)
        f.close()        

        # Log summary
        if self.logNormal():
            self.log("\n")
            self.log.header1("Summary")
            self.log.parformat("Data files found")
            self.log(str(len(files)))
            self.log("\n")
            self.log.parformat("Data files copied")
            self.log(str(n_copied))
            self.log("\n")
            self.log.parformat("Size")
            self.log("%.2f"%(float(total_size)*1e-6)+" MB")
            self.log("\n")
            self.log.parformat("Local configs")
            self.log(str(len(newconfigs)))
            self.log("\n")
            if self.logTerse():
                self.log("\n")
                self.log.header3("Content of master index")
                for config in newconfigs:
                    self.log(str(config["name"]))
                    self.log("\n")
        
        # Return
        return       

    def copy(self, source, clobber):     
        """
        copy file to m_outpath
        """
        # Get file destination
        destination = os.path.join(self.m_outpath,os.path.relpath(source, self.remote_base))
            
        # Initialise return value
        filesize = 0.0
        
        # Logging
        if self.logVerbose():
            self.log("\n")  
            self.log.header3("Copy file")
            self.log.parformat("Source")
            self.log(str(source))
            self.log("\n")
            self.log.parformat("Destination")
            self.log(str(destination))
            self.log("\n")  
            self.log.parformat("Already exists")
            self.log(str(os.path.isfile(destination)))
            self.log("\n")
            self.log.parformat("Overwrite")
            self.log(str(clobber))
            self.log("\n")
            self.log.parformat("Copying")
        
        # Flag if file destination is already available
        is_file = os.path.isfile(destination)
                
        # check if file could be skipped because clobber=no
        if is_file and clobber == False:
            if self.logVerbose():
                self.log("Skip (clobber=no)")
                self.log("\n")
        
        # check if file could be skipped because it is the same file
        elif is_file and os.path.samefile(destination, source):
            if self.logVerbose():
                self.log("Skip (same file)")
                self.log("\n")
            
        else:  
            # Create directory if not existent
            dest_dir = os.path.dirname(destination)         
            if not os.path.isdir(dest_dir):
                os.makedirs(dest_dir)
            
            # Copy file (if source not identical to destination)
            shutil.copy2(source, destination)
            
            # Get filesize
            filesize = os.stat(destination).st_size
            
            # Logging
            if self.logVerbose():
                self.log("Done!")
                self.log("\n")
                   
        # Return 
        return filesize
    
    
    def merge(self, localfits, remotefits, hduname, clobber):
        """
        merge remote and local fits file
        If local fits file not present, a new one is created
        """    
        
        if self.logTerse():
            self.log.parformat("Local file")
            self.log(str(localfits))
            self.log("\n")
            self.log.parformat("Remote file")
            self.log(str(remotefits))
            self.log("\n")
            self.log.parformat("HDU name")
            self.log(str(hduname))
            self.log("\n")
            self.log.parformat("Conflict preference")
            if clobber:
                self.log("Remote")
            else:
                self.log("Local")
            self.log("\n")
        
        # Check if local fits file is available
        if os.path.isfile(localfits):    
            local = gammalib.GFits(str(localfits))
            local_hdu = local[hduname]
        else:
            # Otherwise load remote fits file and delete rows
            local = gammalib.GFits(str(remotefits))
            local_hdu = local[hduname]
            local_hdu.remove_rows(0,local_hdu.nrows())
        
        # find local obs_id 
        local_obs = []
        lobs_col = local_hdu["OBS_ID"]
        for obs_id in lobs_col:
            local_obs.append(obs_id)      
        
        # Create selection expression for opening remote fits file
        selection = ""
        if len(self.runlist):
            selection += "["+hduname+"]["
            for run in self.runlist:
                selection += "OBS_ID=="+str(run)
                selection += "||"
            selection = selection[:-2]
            selection += "]"
        
        # Create file name
        remotefile = remotefits + selection
        
        # open remote fits file
        remote = gammalib.GFits(str(remotefile))
        remote_hdu = remote[hduname]
        
        # Get remote obs_id
        remote_obs = []
        robs_col = remote_hdu["OBS_ID"]
        for obs_id in robs_col:
            remote_obs.append(obs_id)
            
        if self.logTerse():
            self.log.parformat("Remote entries")
            self.log(str(len(remote_obs)))
            self.log("\n")
            self.log.parformat("Local entries")
            self.log(str(len(local_obs)))
            self.log("\n")
        
        # initialise counter for logging
        removed_rows = 0
        
        # If clobber=yes, we use the remote content to overwrite corresponding local content
        if clobber:
            
            # Loop over remote obs_ids
            for remote_obs_id in remote_obs:
                
                # Check if remote obs_id is present locally
                if remote_obs_id in local_obs:
                    
                    # Remove local obs_id entries
                    table_has_obsid = True
                    while table_has_obsid:
                        for i in range(local_hdu.nrows()):
                            if remote_obs_id == local_hdu["OBS_ID"][i]:
                                local_hdu.remove_rows(i,1)
                                removed_rows += 1
                                break
                        table_has_obsid = False
                        for i in range(local_hdu.nrows()):
                            if remote_obs_id == local_hdu["OBS_ID"][i]:
                                table_has_obsid = True
                                break
                            
        # If clobber=false, we keep the local content unchanged and just expand by remote content        
        else:
            # Loop over remote obs_ids
            for local_obs_id in local_obs:
                
                # Check if local obs_id is present remotely
                if local_obs_id in remote_obs:
                    
                    # Remove remote obs_id entries
                    table_has_obsid = True
                    while table_has_obsid:
                        for i in range(remote_hdu.nrows()):
                            if local_obs_id == remote_hdu["OBS_ID"][i]:
                                remote_hdu.remove_rows(i,1)
                                removed_rows += 1
                                break
                        table_has_obsid = False
                        for i in range(remote_hdu.nrows()):
                            if local_obs_id == remote_hdu["OBS_ID"][i]:
                                table_has_obsid = True
                                break
    
        # Store number of local rows
        old_local_rows = local_hdu.nrows()
        
        # Add rmeote HDUs
        local_hdu.insert_rows(old_local_rows, remote_hdu.nrows())
        
        if self.logTerse():
            if clobber:
                self.log.parformat("Removed local rows")
            else:
                self.log.parformat("Skipped remote rows")
            self.log(str(removed_rows))
            self.log("\n")

        # Loop over columns 
        for i in range(local_hdu.ncols()):
            
            # Get local and remote columns
            local_col = local_hdu[i]
            remote_col = remote_hdu[i]    
            
            # Loop over entries and merge
            for j in range(remote_col.length()):
                local_col[j+old_local_rows] = remote_col[j]
        
        # Save local fits file
        local.saveto(str(localfits), True)
               
        # Return
        return
        
# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Copy IACT data
    """
    # Create instance of application
    app = csiactcopy(sys.argv)

    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
