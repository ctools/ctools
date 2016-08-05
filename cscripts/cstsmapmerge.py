#!/usr/bin/env python
# ==========================================================================
# Merge Test Statistic maps
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
import glob
import os
import gammalib
import ctools


# ================== #
# cstsmapmerge class #
# ================== #
class cstsmapmerge(ctools.cscript):
    """
    Merge Test Statistic maps.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name    = 'cstsmapmerge'
        self._version = '1.1.0'

        # Initialise class members
        self._files        = None
        self._in_filename  = ''
        self._tsmap        = gammalib.GSkyMap() # Empty sky map
        self._statusmap    = gammalib.GSkyMap() # Empty sky map
        self._maps         = []
        self._mapnames     = []
        self._merged_files = []
        self._overwrite    = True
        self._delete       = False

        # Initialise application by calling the appropriate class
        # constructor.
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile.
        """
        # Get input maps string
        inmaps = self['inmaps'].string()

        # Handle ASCII files. If the file names given in the ASCII are
        # relative filenames it is assumed that the filename is given
        # relative to the location of the file.
        if '@' == inmaps[0]:
            filename    = inmaps.replace('@','')
            self._files = open(filename).read().splitlines()
            dirname     = os.path.dirname(filename)
            files       = []
            for f in self._files:
                if f[0] != '/':
                    fname = dirname + '/' + f
                else:
                    fname = f
                files.append(fname)
            self._files = files
        
        # Handle wild card strings
        elif '*' in inmaps:
            self._files = glob.glob(inmaps)
        
        # Handle space-separated list
        elif ' ' in inmaps:
            self._files = inmaps.split(' ')
        
        # Handle semi-colon separated list
        elif ';' in inmaps:
            self._files = inmaps.split(';')
            
        # Throw exception if input models cannot be decoded
        else:
            msg = 'Parameter "inmaps" must contain either an @ASCII '\
                  'file, a semi-colon-separated or whitespace-separated '\
                  'list of files or a wildcard string.'
            raise RuntimeError(msg) 

        # Check number of files. We need at least two files.
        if len(self._files) <= 1:
            msg = 'Need at least two files to start merging, '+\
                  str(len(self._files))+' file(s) given.'
            raise RuntimeError(msg) 

        # Get other parameters
        self._overwrite = self['overwrite'].boolean()
        self._delete    = self['delete'].boolean()
        
        # Read ahead output filename
        if self._read_ahead():
            self['outmap'].filename()
        
        # Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return
 
    def _init_ts_map(self, fitsfile):
        """
        Initialise Test Statistic map.
        """
        # Set filename
        self._in_filename = fitsfile
        
        # Open FITS file
        fits = gammalib.GFits(fitsfile)

        # Read TS and status maps
        self._tsmap     = gammalib.GSkyMap()
        self._tsmap.read(fits[0])
        self._statusmap = gammalib.GSkyMap()
        self._statusmap.read(fits['STATUS MAP'])
        
        # Get other maps 
        self._maps     = []
        self._mapnames = []
        
        # Loop over extensions
        for hdu in fits:
            
            # Leave out primary and status extension
            if hdu.extname() != 'IMAGE' and hdu.extname() != 'STATUS MAP':
                
                # Add present maps
                skymap = gammalib.GSkyMap()
                skymap.read(hdu)
                self._maps.append(skymap)
                self._mapnames.append(hdu.extname())
                
        # Close FITS file       
        fits.close()
        
        # Return
        return

    def _merge_ts_map(self, fitsfile):
        """
        Merge TS map from FITS file into output TS map.

        Args:
            fitsfile: FITS file to be merged.
        """
        # Open FITS file
        fits = gammalib.GFits(fitsfile)

        # Read TS and status maps
        add_tsmap     = gammalib.GSkyMap()
        add_tsmap.read(fits[0])
        add_statusmap = gammalib.GSkyMap()
        add_statusmap.read(fits['STATUS MAP'])
        
        # Get other maps 
        add_maps = []
        
        # Loop over extensions
        for hdu in fits:
            
            # Leave out primary and status extension
            if hdu.extname() != 'IMAGE' and hdu.extname() != 'STATUS MAP':
                
                # Add present maps
                skymap = gammalib.GSkyMap()
                skymap.read(hdu)
                add_maps.append(skymap)
                
        # Close FITS file
        fits.close()
        
        # Compare size of maps
        if not len(add_maps) == len(self._maps):
            msg = 'Cannot merge map "'+fitsfile+'" into map "'+\
                  self._in_filename+'" since the number of parameters ' +\
                  'between both maps is different.'
            raise RuntimeError(msg) 

        # Loop over bins    
        for i in range(self._tsmap.npix()):

            # Consider only bins that have been computed
            if add_statusmap[i] > -0.5:
                
                # Raise exception if this bin has already been computed
                if self._statusmap[i] > -0.5 and not self._overwrite:
                    msg = 'Attempt to merge bin which apparently has '+\
                          'already been merged. File "'+fitsfile+'" '+\
                          'contains already merged bins. Set hidden '+\
                          'parameter "overwrite=yes" to avoid this error.'
                    raise RuntimeError(msg)
                
                # Copy TS values
                self._tsmap[i] = add_tsmap[i]
                
                # Copy status 
                self._statusmap[i] = add_statusmap[i]
                
                # Loop over maps and copy entries
                for j in range(len(self._maps)):
                    self._maps[j][i] = add_maps[j][i]
        
        # Return
        return

    def _get_number_of_ts_pixels(self, fitsfile):
        """
        Return number of pixels with TS values.

        Args:
            fitsfile: FITS file to be merged.

        Returns:
            Number of pixels for which TS has been computed.
        """
        # Get status map for this file
        status = gammalib.GSkyMap(fitsfile+'[STATUS MAP]')

        # Count number of pixels with TS status set
        count = 0
        for pix in status:
            if pix > -0.5:
                count += 1

        # Return number of pixels
        return count


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
        
        # Write header into logger
        self._log_header1(gammalib.TERSE, 'Merge TS maps')
        
        # Initialise file to start with
        file0              = ''
        self._merged_files = []
        
        # Test files for the entry status map. Use the first one to appear
        # useful
        for fitsfile in self._files:

            # Skip file if it's not a FITS file
            if not gammalib.GFilename(fitsfile).is_fits():
                
                # Log message
                self._log_value(gammalib.EXPLICIT, 'Skip file', fitsfile +' (not a FITS file)')

                # Continue
                continue

            # Open FITS file
            fits = gammalib.GFits(fitsfile)

            # If file contains a status map then use it
            if fits.contains('STATUS MAP'):
                
                # Close FITS file
                fits.close()
                
                # Set filename to initial filename
                file0 = fitsfile
                
                # Log message
                pix_info = ' (%d TS pixels computed)' % self._get_number_of_ts_pixels(fitsfile)
                self._log_value(gammalib.TERSE, 'Initial TS map file', fitsfile + pix_info)
                
                # Leave loop
                break

            # ... otherwise signal that file is useless
            else:
                
                # Close fits file
                fits.close()  
                
                # Info message             
                self._log_value(gammalib.EXPLICIT, 'Skip file', fitsfile +' (no "STATUS MAP" extension)')
                
                # Continue
                continue

        # Signal if no suitable file was found
        if file0 == '':
            msg = 'None of the provided files seems to be a sliced ' + \
                          'TS map file (none has a "STATUS MAP" ' + \
                          'extension)'
            raise RuntimeError(msg)
                          

        # ... otherwise merge files
        else:
            
            # Copy file list
            workfiles = self._files
            
            # Remove entry which will be used to initalise the map
            workfiles.remove(file0)
            
            # Initialise map from first file
            self._init_ts_map(file0) 
            
            # Append to added files 
            self._merged_files.append(file0)
            
            # Loop over files
            for fitsfile in workfiles:
                
                # Skip if file is not FITS
                if not gammalib.GFilename(fitsfile).is_fits():
                    self._log_value(gammalib.EXPLICIT, 'Skip file', fitsfile +' (not a FITS file)')
                    continue   

                # Open FITS file
                fits = gammalib.GFits(fitsfile)
                
                # Skip if file does not contain status map
                if not fits.contains('STATUS MAP'):
                    fits.close()
                    self._log_value(gammalib.EXPLICIT, 'Skip file', fitsfile +' (no "STATUS MAP" extension)')
                    continue

                # Close FITS file
                fits.close()
                
                # Logging
                pix_info = ' (%d TS pixels computed)' % self._get_number_of_ts_pixels(fitsfile)
                self._log_value(gammalib.TERSE, 'Merge TS map file', fitsfile + pix_info)
                    
                # Merge TS map
                self._merge_ts_map(fitsfile)
                
                # Append FITS file to merged files
                self._merged_files.append(fitsfile)
                
        # Return
        return

    def save(self):
        """ 
        Save TS map and remove slices if requested.
        """

        # Get output filename in case it was not read ahead
        outmap = self['outmap'].filename()
        
        # Create FITS file
        fits = gammalib.GFits()
        
        # Write TS map into primary
        self._tsmap.write(fits)

        # Loop over maps and write them to fits
        for i in range(len(self._maps)):
            self._maps[i].write(fits)
        
        # Set map names as extensions
        for i in range(len(self._mapnames)):   
            fits[i+1].extname(self._mapnames[i])

        # Treat status map separately 
        self._statusmap.write(fits)
        fits[fits.size()-1].extname('STATUS MAP')
        
        # Check if map is fully done
        bins_merged = 0
        for pix in self._statusmap:
            if pix > -0.5:
                bins_merged += 1

        # Set boolean if map is fully done
        done = (bins_merged == self._tsmap.npix())
        
        # Log summary header
        self._log_header2(gammalib.TERSE, 'Merging Summary')
                
        # Write computing status if we are not done yet
        if not done:
            
            # Log bins that were computed and those that were missing
            self._log_value(gammalib.TERSE, 'TS map bins', self._tsmap.npix())
            self._log_value(gammalib.TERSE, 'Bins computed and merged', bins_merged )
            self._log_value(gammalib.TERSE, 'Bins missing', self._tsmap.npix() - bins_merged)
        
        else:
            
            # Write success message
            self._log_string(gammalib.TERSE, 'TS map was fully computed and successfully merged') 
        
        # Write header
        self._log_header1(gammalib.TERSE, 'Save TS map')
        
        # Log filename
        self._log_value(gammalib.TERSE, 'TS map file', outmap.url())
        
        # Save FITS file
        fits.saveto(outmap, self._clobber()) 

        # Delete TS input maps if requested and map is fully done
        if self._delete and done:
            
            # Log header
            self._log_header2(gammalib.TERSE, 'Delete slices TS map files')
            
            # Loop over sliced files
            for filename in self._merged_files:
                
                # Remove file
                os.remove(filename)
                
                # Log message about file
                self._log_value(gammalib.TERSE,'Deleted input file', filename)
        
        # Return
        return

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

        # Save TS map if required
        self.save()

        # Return
        return    
        

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cstsmapmerge(sys.argv)
    
    # Execute application
    app.execute()
