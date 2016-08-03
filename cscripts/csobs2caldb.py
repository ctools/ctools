#!/usr/bin/env python
# ==========================================================================
# Creates a calibration database entry for an IACT observation.
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
import os
import datetime
import gammalib
import ctools
from cscripts import calutils


# ================= #
# csobs2caldb class #
# ================= #
class csobs2caldb(ctools.cscript):
    """
    Creates a calibration database entry for an IACT observation.
    
    The creation of a calibration database entry is useful for performing
    simulations for current Imaging Air Cherenkov Telescopes (IACTs).
    The class takes an observation definition XML file as input and uses
    one observation to create a calibration database entry. By default
    the first observation will be used, but the user can specify the
    index of any observation using the hidden "index" parameter.
    """
    
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name and version
        self._name    = 'csobs2caldb'
        self._version = '1.2.0'

        # Initialise members
        self._observation = gammalib.GCTAObservation()
        self._mission     = 'cta'
        self._caldb       = 'cta'
        self._outfile     = gammalib.GFilename('irf_file.fits')
        self._base_dir    = ''
        self._cal_dir     = ''
        self._rsp_dir     = ''
        self._caldb_inx   = gammalib.GFits()
        self._irf_fits    = gammalib.GFits()

        # Initialise observation container from constructor arguments.
        self._obs, argv = self._set_input_obs(argv)
        
        # Initialise script by calling the appropriate class constructor.
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile

        Raises:
            ValueError, IndexError, & RuntimeError.
        """
        # Load observation container from "inobs" parameter in case that
        # the observation container is still empty.
        if self._obs.is_empty():
            self._obs = gammalib.GObservations(self['inobs'].filename())

        # Raise an exception if the observation container is empty
        if self._obs.is_empty():
            raise ValueError('No or empty observation provided for '
                             'parameter "inobs".')
        
        # Get index of observation to be used
        index = self['index'].integer()

        # Raise an exception if the index is not valid
        if index < 0 or index >= self._obs.size():
            raise IndexError('Parameter "index=%d" outside the validity '
                             'range [0,%d].' % (index, self._obs.size()))
        
        # Get observation
        self._observation = self._obs[index]

        # Make sure we have a CTA observation
        if not self._observation.classname() == 'GCTAObservation':
            raise RuntimeError('Input observation not of type '
                               '"GCTAObservation".')
        
        # Make sure we have an IRF response associated with the CTA
        # observation
        if not self._observation.response().classname() == 'GCTAResponseIrf':
            raise RuntimeError('Response of input observation not of '
                               'type "GCTAResponseIrf".')
        
        # Get calibration database name. If the "caldb" parameter is "NONE"
        # or empty then use the instrument name from the observation as
        # calibration database name.
        self._caldb = self['caldb'].string()
        if (self._caldb == 'NONE' or len(self._caldb) == 0):
            self._caldb = self._observation.instrument().lower()

        # Get instrument response function file name
        self._outfile = self['outfile'].filename()

        # Make sure that remaining user parameters are queried now. We
        # do not store the actual parameter values as we do not want
        # too many instance attributes with enhances the maintenance
        # costs.
        self['irf'].string()
        self['rootdir'].string()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _make_irf_file(self):
        """
        Creates an IRF FITS file
        """
        # Write header into logger
        if self._logTerse():
            self._log("\n")
            self._log.header2("Creating IRF file")

        # Get response for the observation
        rsp = self._observation.response()
        
        # Extract response file names
        fname_aeff  = rsp.aeff().filename()
        fname_psf   = rsp.psf().filename()
        fname_edisp = rsp.edisp().filename()
        fname_bkg   = rsp.background().filename()

        # Log filenames
        self._log_header3(gammalib.NORMAL, 'IRF input files')
        self._log_value(gammalib.NORMAL, 'Effective area', fname_aeff.url())
        self._log_value(gammalib.NORMAL, 'Point spread function', fname_psf.url())
        self._log_value(gammalib.NORMAL, 'Energy dispersion', fname_edisp.url())
        self._log_value(gammalib.NORMAL, 'Background rate', fname_bkg.url())
    
        # Open FITS files of response components
        fits_aeff  = gammalib.GFits(fname_aeff)
        fits_psf   = gammalib.GFits(fname_psf)
        fits_edisp = gammalib.GFits(fname_edisp)
        fits_bkg   = gammalib.GFits(fname_bkg)

        # Get extension names
        ext_aeff  = fname_aeff.extname("EFFECTIVE AREA")
        ext_psf   = fname_psf.extname("POINT SPREAD FUNCTION")
        ext_edisp = fname_edisp.extname("ENERGY DISPERSION")
        ext_bkg   = fname_bkg.extname("BACKGROUND")
        
        # Create empty FITS file
        fits = gammalib.GFits()
        
        # Append IRF component to FITS file
        fits.append(fits_aeff[ext_aeff])
        fits.append(fits_psf[ext_psf])
        fits.append(fits_edisp[ext_edisp])
        fits.append(fits_bkg[ext_bkg])

        # Log resulting FITS file
        if self._logNormal():
            self._log(str(fits))
            self._log("\n")
        if self._logExplicit():
            self._log(str(fits[ext_aeff].header()))
            self._log("\n")
            self._log(str(fits[ext_psf].header()))
            self._log("\n")
            self._log(str(fits[ext_edisp].header()))
            self._log("\n")
            self._log(str(fits[ext_bkg].header()))
            self._log("\n")
        
        # Return fits file
        return fits

    def _make_caldb_index(self):
        """
        Creates an IRF index FITS file
        """
        # Write header into logger
        if self._logTerse():
            self._log('\n')
            self._log.header2('Creating database index file')

        # Open calibration database index
        indx_file = self._base_dir + '/caldb.indx'
        
        # Open index file (or create one if it does not exist)
        cif = gammalib.GFits(indx_file, True)
        
        # Retrieve "CIF" table
        if cif.contains('CIF'):
            table = cif.table('CIF')

        # ... or create binary table if no "CIF" table exists yet, append
        # binary table to CIF file and return a reference to the appended
        # table to work with
        else:       
            table = cif.append(calutils.create_cif_table())

        # Check if output config already exist
        has_config = False
        row_index = -1
        for caldir in table['CAL_DIR']:
            row_index += 1
            if caldir == self._cal_dir:
                has_config = True
                break
         
        # Create columns if not available   
        if not has_config:           
            # Append 4 rows to CIF extension
            table.append_rows(4)
            row_index = table.nrows()
        else:
            row_index += 4

        # Add generic information for these 4 rows
        for i in range(4):

            # Set row number
            row = i + row_index-4

            # Set date
            now = str(datetime.datetime.now())
            today, time = now.split()
            
            # Set element
            table['TELESCOP'][row] = self._mission
            table['INSTRUME'][row] = self._caldb
            table['DETNAM'][row]   = 'NONE'
            table['FILTER'][row]   = 'NONE'
            table['CAL_DEV'][row]  = 'ONLINE'
            table['CAL_CLAS'][row] = 'BCF'
            table['CAL_DTYP'][row] = 'DATA'
            table['CAL_VSD'][row]  = today
            table['CAL_VST'][row]  = time.split('.')[0]
            table['REF_TIME'][row] = 51544.0
            table['CAL_QUAL'][row] = 0
            table['CAL_CBD'][row]  = 'NAME('+self['irf'].string()+')'
            table['CAL_DATE'][row] = today.replace('-','/')[2:]
            table['CAL_DIR'][row]  = self._cal_dir
            table['CAL_FILE'][row] = self._outfile

        # Add effective area information
        row = row_index-4
        table['CAL_CNAM'][row] = 'EFF_AREA'
        table['CAL_DESC'][row] = self._caldb+' effective area'

        # Add point spread function information
        row = row_index-3
        table['CAL_CNAM'][row] = 'RPSF'
        table['CAL_DESC'][row] = self._caldb+' point spread function'
        
        # Add energy dispersion information
        row = row_index-2
        table['CAL_CNAM'][row] = 'EDISP'
        table['CAL_DESC'][row] = self._caldb+' energy dispersion'

        # Add background information
        row = row_index-1
        table['CAL_CNAM'][row] = 'BGD'
        table['CAL_DESC'][row] = self._caldb+' background'

        # Log resulting FITS table and header
        if self._logNormal():
            self._log(str(table))
            self._log('\n')
        if self._logExplicit():
            self._log(str(table.header()))
            self._log('\n')
        
        # Return CIF FITS file
        return cif
        
    def _make_dirs(self):
        """
        Make CALDB directories
        """
        # Write header into logger
        if self._logTerse():
            self._log("\n")
            self._log.header2("Creating directory structure")

        # Create calibration database        
        caldb = gammalib.GCaldb(self["rootdir"].string())
        
        # Set calibration directory 
        self._cal_dir  = "data"
        self._cal_dir += "/"+self._mission.lower()
        self._cal_dir += "/"+self._caldb.lower()
        self._cal_dir += "/bcf/"+self["irf"].string()
        
        # Set absolute path
        self._base_dir = caldb.rootdir() +"/data"
        self._base_dir += "/"+self._mission.lower()
        self._base_dir += "/"+self._caldb.lower()
        
        # Set directory for irf file
        self._rsp_dir = caldb.rootdir() + "/" + self._cal_dir

        # Log resulting FITS table
        self._log_value(gammalib.NORMAL, 'Calibration directory', self._cal_dir)
        self._log_value(gammalib.NORMAL, 'Base directory', self._base_dir)
        if not os.path.isdir(self._rsp_dir):
            name = 'IRF directory'
        else:
            name = 'IRF directory (existing)'
        self._log_value(gammalib.NORMAL, name, self._rsp_dir)
        
        # Create IRF directory is it does not yet exist
        if not os.path.isdir(self._rsp_dir):
            os.makedirs(self._rsp_dir)

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
        
        # Write input parameters and header into logger
        if self._logTerse():
            self._log('\n')
            self._log.header1('Creating CALDB entry')
        
        # Create directory structure
        self._make_dirs()
         
        # Create/update calibration database   
        self._caldb_inx = self._make_caldb_index()
        
        # Create response file
        self._irf_fits = self._make_irf_file()
                
        # Return
        return

    def save(self):
        """
        Save calibration database FITS file
        """
        # Write header into logger
        self._log_header1(gammalib.TERSE, 'Save calibration database')

        # Set response filename
        filename = self._rsp_dir + '/' + self._outfile
        
        # Write filenames into logger
        self._log_value(gammalib.NORMAL, 'CALDB index file',
                        self._caldb_inx.filename().url())
        self._log_value(gammalib.NORMAL, 'Response file', filename)

        # Save caldb index file
        self._caldb_inx.save(self._clobber())

        # Save response file
        self._irf_fits.saveto(filename, self._clobber())
        
        # Return
        return
        
    def execute(self):
        """
        Execute the script
        """
        # Open logfile
        self.logFileOpen()

        # Run the script
        self.run()

        # Save to disk
        self.save()

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csobs2caldb(sys.argv)
    
    # Execute application
    app.execute()
