#!/usr/bin/env python
# ==========================================================================
# Generates an IACT observation definition XML file.
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
import json
import gammalib
import ctools


# =============== #
# csiactobs class #
# =============== #
class csiactobs(ctools.cscript):
    """
    Generates an IACT observation definition XML file
    
    This class implements the creation of a observation xml file for IACT
    data analysis. This class is dedicated for use inside a IACT
    Collaboration, i.e. it can only be used if you have access to IACT data
    in FITS format. The FITS data has to be structured and in the format
    described here:
    http://gamma-astro-data-formats.readthedocs.org/en/latest/
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name and version
        self._name    = 'csiactobs'
        self._version = '1.2.0'

        # Initialise some members
        self._obs              = gammalib.GObservations()
        self._ebounds          = gammalib.GEbounds()
        self._datapath         = os.getenv('VHEFITS','')
        self._prodname         = ''
        self._xml              = gammalib.GXml()
        self._models           = gammalib.GModels()
        self._runlist          = []
        self._runlistfile      = gammalib.GFilename()
        self._bkgpars          = 0
        self._master_indx      = ''
        self._use_bkg_scale    = False
        self._ev_hiera         = ['']
        self._aeff_hiera       = ['']
        self._psf_hiera        = ['']
        self._bkg_hiera        = ['']
        self._edisp_hiera      = ['']
        self._bkg_mod_hiera    = ['']
        self._bkg_gauss_norm   = 1.0
        self._bkg_gauss_index  = 0.0
        self._bkg_gauss_sigma  = 1.0
        self._bkg_aeff_index   = 0.0
        self._bkg_aeff_norm    = 1.0
        self._bkg_range_factor = 1.0
        self._hdu_index        = ''
        self._obs_index        = ''
        self._subdir           = ''
        self._debug            = False

        # Initialise empty observation definition XML file
        self._xml.append(gammalib.GXmlElement('observation_list '
                                              'title="observation list"'))

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Append an observation list to XML instance
        self._xml.append(gammalib.GXmlElement('observation_list title="observation list"'))
        self._xml_obslist = self._xml.element('observation_list', 0)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation
        """
        # Set data path
        if self._datapath == '':
            self._datapath = self['datapath'].string()
        
        # Query input parameters
        self['inmodel'].filename()

        # Expand environment
        self._datapath = gammalib.expand_env(self._datapath)
        
        # Read FITS production
        self._prodname = self['prodname'].string()
        
        # Read runlist file if list not already filled
        if len(self._runlist) == 0:
            
            # Get file name
            self._runlistfile = self['infile'].filename()
            
            # Read runlist from file
            runfile = open(self._runlistfile.url())
            for line in runfile.readlines():
                if len(line) == 0:
                    continue
                if line[0] == '#':
                    continue
                if len(line.split()) > 0:
                    self._runlist.append(line.split()[0])
            runfile.close()
        
        # Read number of background parameters
        self._bkgpars = self['bkgpars'].integer()

        # Query ahead output parameters
        if self._read_ahead():
            self['outmodel'].filename()
            self['outobs'].filename()
        
        # Master index file name
        self._master_indx = self['master_indx'].string()
        
        # Read flag for background scaling factor
        self._use_bkg_scale = self['bkg_scale'].boolean()
        
        # Read hierarchy of file loading
        self._ev_hiera      = self['ev_hiera'].string().split('|')
        self._aeff_hiera    = self['aeff_hiera'].string().split('|')
        self._psf_hiera     = self['psf_hiera'].string().split('|')
        self._bkg_hiera     = self['bkg_hiera'].string().split('|')
        self._edisp_hiera   = self['edisp_hiera'].string().split('|')
        self._bkg_mod_hiera = self['bkg_mod_hiera'].string().split('|')
        
        # Read hidden background parameters
        self._bkg_gauss_norm   = self['bkg_gauss_norm'].real()
        self._bkg_gauss_index  = self['bkg_gauss_index'].real()
        self._bkg_gauss_sigma  = self['bkg_gauss_sigma'].real()
        self._bkg_aeff_index   = self['bkg_aeff_index'].real()
        self._bkg_aeff_norm    = self['bkg_aeff_norm'].real()
        self._bkg_range_factor = self['bkg_range_factor'].real()
        
        # Open master index file and look for prodname 
        master_file = os.path.join(self._datapath, self._master_indx)
        if not os.path.isfile(master_file):
            raise RuntimeError('FITS data store not available. No master index file found. Make sure the file is copied from the server and your datapath is set correctly.')

        # Open and load JSON file
        json_data = open(master_file).read()
        data      = json.loads(json_data)    
        if not 'datasets' in data:
            raise RuntimeError('Key "datasets" not available in master index file.')

        # Get configurations
        configs = data['datasets']

        # Initialise HDUs
        self._hdu_index = self._obs_index = ''

        # Get HDUs
        for config in configs:
            
            # Check if prodname is present
            if self._prodname == config['name']:
                self._hdu_index = str(os.path.join(self._datapath,
                                                   config['hduindx']))
                self._obs_index = str(os.path.join(self._datapath,
                                                   config['obsindx']))
                
                # Leave loop if index file names were found
                break

        # Check index files
        if self._hdu_index == '' or self._obs_index == '':
            raise RuntimeError('*** ERROR: FITS data store "'+self._prodname+'" not available. Run csiactdata to get a list of available storage names')
        
        # Check HDU names
        filename = gammalib.GFilename(self._hdu_index+'[HDU_INDEX]')
        if not filename.is_fits():
            raise RuntimeError('*** ERROR: HDU index file "'+self._hdu_index+'[HDU_INDEX]" for FITS data store "'+self._prodname+'" not available. Check your master index file or run csiactdata to get a list of available storage names.')

        # Check for existence of 'BKG_SCALE' in the observation index file if required
        if self._use_bkg_scale:
            
            # Create filename
            filename = gammalib.GFilename(self._obs_index+'[OBS_INDEX]')
            
            # Check if it is a FITS file
            if filename.is_fits():
                
                # Open FITS file
                fits = gammalib.GFits(self._obs_index)
                
                # Check if column "BKG_SCALE" is found and signal its possible usage
                if not fits['OBS_INDEX'].contains('BKG_SCALE'):
                    self._use_bkg_scale = False
                    
                # Close FITS file
                fits.close()
                
            else:
                # Signal that there is no background scale
                self._use_bkg_scale = False
        
        # Create base data directory from hdu index file location
        self._subdir  = os.path.dirname(self._hdu_index)  
        self._debug   = False # Debugging in client tools

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return
    
    def _background_spectrum(self, prefactor, index, emin = 0.01, emax = 100.0):
        """
        Create a background spectrum model dependent on user parameters
        
        Parameters
        ----------
        prefactor : float
            Power law prefactor of spectral model
        index : float
            Power law index of spectral model
        emin : float
            Minimum energy (in case a spectral node function is required)
        emax : float
            Maximum energy (in case a spectral node function is required)
        
        Returns
        -------
        spec : `~gammalib.GModelSpectral()`
            Spectral model for the background shape     
        """
        # Handle constant spectral model 
        if index == 0.0 and self._bkgpars <= 1:
            spec = gammalib.GModelSpectralConst()
            
            # Set parameter range
            spec['Normalization'].min(prefactor / self._bkg_range_factor)
            spec['Normalization'].max(prefactor * self._bkg_range_factor)
            spec['Normalization'].value(prefactor)
            
            # Freeze or release normalisation parameter
            if self._bkgpars == 0:
                spec['Normalization'].fix()
            else:
                spec['Normalization'].free()
                
        else:
                
            # Create power law model 
            if self._bkgpars <= 2:
                
                # Set Power Law model with arbitrary pivot energy
                pivot = gammalib.GEnergy(1.0,'TeV')
                spec  = gammalib.GModelSpectralPlaw(prefactor, index, pivot)  
                 
                # Set parameter ranges
                spec[0].min(prefactor / self._bkg_range_factor)
                spec[0].max(prefactor * self._bkg_range_factor)
                spec[1].scale(1)
                spec[1].min(-5.0)
                spec[1].max(5.0)
                
                # Set number of free parameters
                if self._bkgpars == 0:
                    spec[0].fix()
                    spec[1].fix()
                elif self._bkgpars == 1:
                    spec[0].free()
                    spec[1].fix()
                else:
                    spec[0].free()
                    spec[1].free()
            
            else:
                
                # Create reference powerlaw
                pivot = gammalib.GEnergy(1.0,'TeV')
                plaw = gammalib.GModelSpectralPlaw(prefactor, index, pivot) 
                
                # Create spectral model and energy values
                spec = gammalib.GModelSpectralNodes()
                
                # Create logarithmic energy nodes
                bounds = gammalib.GEbounds(self._bkgpars,
                                           gammalib.GEnergy(emin,'TeV'),
                                           gammalib.GEnergy(emax,'TeV'), True)
                
                # Loop over bounds and set intensity value
                for i in range(bounds.size()):     
                    energy = bounds.elogmean(i)
                    value = plaw.eval(energy)
                    
                    # Append energy, value - tuple to node function
                    spec.append(energy, value)
                
                # Loop over parameters
                for par in spec:
                    
                    # Fix energy nodes
                    if 'Energy' in par.name():
                        par.fix()
                    
                    # Release intensity nodes
                    elif 'Intensity' in par.name():   
                        value = par.value() 
                        par.scale(value)    
                        par.min(value / self._bkg_range_factor)
                        par.max(value * self._bkg_range_factor)
                        
        # Return spectrum
        return spec

    def _iact_background(self, telescope, obs_id, bkg_scale, bkgtype,
                         emin=0.01, emax=100):
        """
        Create an IACT background model
        
        Parameters
        ----------
        telescope : str
            Name of telescope
        obs_id : str
            Observation ID
        bkg_scale : float
            Background scaling factor
        bkgtype : str
            Type of background (irf,aeff, or gauss)
        
        Returns
        -------
        model : `~gammalib.GModelData()`
            Background model for IACT observation
        """
        # Handle IrfBackground
        if bkgtype == 'irf':
            
            # Set parameters to have a constant spectral model
            prefactor  = 1.0
            index      = 0.0
            prefactor *= bkg_scale
            
            # Create background spectrum
            spec = self._background_spectrum(prefactor, index, emin, emax)
            
            # Create background model
            bck = gammalib.GCTAModelIrfBackground(spec)
        
        # Set AeffBackground   
        elif bkgtype == 'aeff':

            # Set background components
            prefactor = bkg_scale * self._bkg_aeff_norm
            spec      = self._background_spectrum(prefactor,
                                                  self._bkg_aeff_index,
                                                  emin, emax)
                
            # Create background model
            bck = gammalib.GCTAModelAeffBackground(spec)
            
        # Set Gaussian Background
        elif bkgtype == 'gauss':

            # Set background components
            prefactor = bkg_scale * self._bkg_gauss_norm
            spec      = self._background_spectrum(prefactor,
                                                  self._bkg_gauss_index,
                                                  emin, emax)
            radial    = gammalib.GCTAModelRadialGauss(self._bkg_gauss_sigma)

            # Create background model
            bck = gammalib.GCTAModelRadialAcceptance(radial, spec)
        
        else:
            msg = 'Background type "'+bkgtype+'" unsupported'
            raise RuntimeError(msg)
        
        # Copy model
        model = bck.clone()
            
        # Assign specific run id
        model.ids(str(obs_id))
        
        # Assign instrument
        model.instruments(telescope)
        
        # Set name (arbitrary)
        model.name('bkg_'+str(obs_id))  
        
        # Turn off TS calculation for background model
        model.tscalc(False)
        
        # Return model
        return model

    def _append_inmodels(self):
        """
        Append input models
        """
        # If there are models provided by "inmodels" then append them now to
        # the model container
        filename = self['inmodel'].filename().url()
        if filename != '' and filename != 'NONE':

            # Log header and input model file
            self._log_header1(gammalib.TERSE, 'Append input models')
            self._log_value(gammalib.NORMAL, 'Input model file', filename)

            # Load input models
            models = gammalib.GModels(filename)

            # Loop over input models and append them to the model container
            for model in models:
                self._models.append(model)
                self._log_value(gammalib.NORMAL, 'Append model',
                                                 '"'+model.name()+'"')

        # Return
        return

    def _write_summary(self):
        """
        Write observation summary
        """
        # Log header
        self._log_header1(gammalib.NORMAL, 'Observation summary')

        # Set energy range dependent on whether boundaries exist or
        # not
        if self._ebounds.size() > 0:
            erange = '%s - %s' % (str(self._ebounds.emin()),
                                  str(self._ebounds.emax()))
        else:
            erange = 'not available'
        
        # Log energy range
        self._log_value(gammalib.NORMAL, 'Energy range', erange)

        # Return
        return

    def _get_filename(self, hdu, hierarchy, formats): 
        """
        Retrieves a filename from a hdu index
        
        Parameters
        ----------
        hdu : `~gammalib.GFitsTable`
            HDU FITS table
        hierarchy : list of str
            List of strings containing the hierarchy how to look for files
        formats : list of str
            File formats available for this observation
        
        Returns
        -------
        filename : str
            Filename of requested file
        """
        
        # Initialise filename
        filename = ''
        
        # Handle hierarchies
        index = -1
        
        # Loop over possible formats
        for version in hierarchy:
            
            # Count number of same available  versions
            n = formats.count(version) 
            
            # Leave loop and store index of first occurence
            if n > 0:
                index = formats.index(version)
                break
        
        # Build filename from information stored in hduindx_hdu
        if index >= 0:
            filename  = os.path.join(self._subdir,
                                     hdu['FILE_DIR'][index],
                                     hdu['FILE_NAME'][index])
            filehdu   = hdu['HDU_NAME'][index]
            filename += '['+filehdu+']'

        # Return filename
        return filename
    
    
    def _is_present(self, filename, obs_id, filetype):
        """
        Checks if a filename is present
        
        Parameters
        ----------
        filename : str
            Filename
        obs_id : int
            Observation ID
        filetype : str
            Type of file
        
        Returns
        -------
        present : bool
            True if filename is present, False otherwise
        """
                
        # Initialise return value
        present = True
        
        # Create filename instance
        fname = gammalib.GFilename(filename)
        
        # Check if file is empty
        if fname.is_empty():
            msg = 'Skipping observation "%s": No %s found' % (obs_id, filetype)
            self._log_string(gammalib.NORMAL, msg)
            present = False
        
        # Check if file exists
        elif not fname.exists():
            msg = ('Skipping observation "%s": %s "%s" does not exist' %
                   (obs_id, filetype, filename))
            self._log_string(gammalib.NORMAL, msg)
            present = False
        
        # Return present flag
        return present    


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

        # Clear models
        self._models.clear()
        
        # Log output
        msg = 'Looping over %d runs' % len(self._runlist)
        self._log_header1(gammalib.TERSE, msg)

        # Initialise energy range values for logging
        self._ebounds.clear()

        # Loop over runs
        for obs_id in self._runlist:

            # Create selection string
            obs_selection = '[OBS_ID=='+str(obs_id)+']'
            
            # Open HDU index file
            hduindx     = gammalib.GFits(self._hdu_index+'[HDU_INDEX]'+
                                         obs_selection)
            hduindx_hdu = hduindx['HDU_INDEX']
            
            # Initialise types and formats
            formats = []
            for i in range(hduindx_hdu.nrows()):  
                formats.append(hduindx_hdu['HDU_CLASS'][i])
            
            # Retrieve filenames            
            eventfile = self._get_filename(hduindx_hdu, self._ev_hiera, formats)
            aefffile  = self._get_filename(hduindx_hdu, self._aeff_hiera, formats)
            psffile   = self._get_filename(hduindx_hdu, self._psf_hiera, formats)
            edispfile = self._get_filename(hduindx_hdu, self._edisp_hiera, formats)
            bkgfile   = self._get_filename(hduindx_hdu, self._bkg_hiera, formats)
            
            # Check for presence of required event file
            if not self._is_present(eventfile, obs_id, "event file"):
                continue
            
            # Check for presence of required effective area file
            if not self._is_present(aefffile, obs_id, "effective area file"):
                continue
            
            # Check for presence of required psf file
            if not self._is_present(psffile, obs_id, "PSF file"):
                continue
            
            # Check for optional energy dispersion file
            if not gammalib.GFilename(edispfile).exists():
                
                # Print warning that edisp cannot be used for this observation
                msg = ('Warning: observation "%s" has no energy dispersion '
                       'information' % obs_id)
                self._log_string(gammalib.NORMAL, msg)
                
                # Set energy dispersion file as empty
                edispfile = ''
            
            # Create a copy of background model hierarchy for this run
            run_bkg_mod_hierarchy = list(self._bkg_mod_hiera)
            
            # Check if background file is available
            # Remove IRF background from hierarchy if file doesnt exist
            if not gammalib.GFilename(bkgfile).exists():
                
                # Set background file as empty
                bkgfile = ''
                
                # Check for IRF background in background model hierarchy
                if 'irf' in run_bkg_mod_hierarchy:
                    
                    # Remove IRF background if file doesnt exist
                    run_bkg_mod_hierarchy.remove('irf')
                
                    # Print warning that edisp cannot be used for this
                    # observation
                    msg = ('Warning: observation "%s" has no background '
                           'information (IRF background cannot be used)' %
                           obs_id)
                    self._log_string(gammalib.NORMAL, msg)
                
                # Check if background model hierarchy is empty
                if len(run_bkg_mod_hierarchy) == 0:
                    
                    # Skip observation if no background can be used
                    msg = ('Skipping observation "%s": No background can be '
                           'used' % obs_id)
                    self._log_string(gammalib.NORMAL, msg)
                    continue
                
                else:
                    # Log if we fall back to next background approach
                    msg = ('Observation "%s": Falling back to background "%s"' %
                           (obs_id, run_bkg_mod_hierarchy[0]))
                    self._log_string(gammalib.NORMAL, msg)

            # Close hdu index file
            hduindx.close()
            
            # Initialise background scale
            bkg_scale = 1.0
            
            # Check if background scale should be used
            if self._use_bkg_scale:
                
                # Read background scale from fits file
                obsindx   = gammalib.GFits(self._obs_index+'[OBS_INDEX]'+
                                           obs_selection)
                bkg_scale = obsindx['OBS_INDEX']['BKG_SCALE'][0]
                obsindx.close()
            
            # Open event FITS file
            fits        = gammalib.GFits(eventfile)
            
            # Get object and telescope strings from event hdu
            events      = fits[gammalib.GFilename(eventfile).extname()]
            object_name = events.string('OBJECT')
            telescope   = events.string('TELESCOP')

            # Close FITS file
            fits.close()
            
            # Open effective area file to look for threshold
            aeff_fits  = gammalib.GFits(aefffile)
            aeff_table = aeff_fits[gammalib.GFilename(aefffile).extname()]
            
            # Set energy range from header keyword if present
            if aeff_table.has_card('LO_THRES') and aeff_table.has_card('HI_THRES'):
                
                # Get safe energy range
                run_emin = gammalib.GEnergy(aeff_table.real('LO_THRES'), 'TeV')
                run_emax = gammalib.GEnergy(aeff_table.real('HI_THRES'), 'TeV')
                
                # Append to ebounds
                self._ebounds.append(run_emin, run_emax)
            
            else:
                # Set default values for energy range
                run_emin = gammalib.GEnergy(10,'GeV')
                run_emax = gammalib.GEnergy(100,'TeV')
            
            # Close Aeff fits file
            aeff_fits.close()     
            
            # Append instrumental background model
            self._models.append(self._iact_background(telescope,
                                                      obs_id,
                                                      bkg_scale,
                                                      run_bkg_mod_hierarchy[0],
                                                      run_emin.TeV(),
                                                      run_emax.TeV()))

            # Append observation to XML and set attributes
            obs = self._xml_obslist.append('observation')
            obs.attribute('name', object_name)
            obs.attribute('id', str(obs_id))
            obs.attribute('instrument', telescope)

            # Append event file
            ev = gammalib.GXmlElement('parameter name="EventList"')
            ev.attribute('file', eventfile)

            # Append effective area
            aeff = gammalib.GXmlElement('parameter name="EffectiveArea"')
            aeff.attribute('file', aefffile)

            # Append PSF
            psf = gammalib.GXmlElement('parameter name="PointSpreadFunction"')
            psf.attribute('file', psffile)

            # Append energy dispersion
            edisp = gammalib.GXmlElement('parameter name="EnergyDispersion"')
            edisp.attribute('file', edispfile)

            # Append background
            bck = gammalib.GXmlElement('parameter name="Background"')
            bck.attribute('file', bkgfile)

            # Assign events and IRFs to observation
            obs.append(ev)
            obs.append(aeff)
            obs.append(psf)
            obs.append(edisp)
            obs.append(bck)

            # Log the observation ID and object name that has been appended
            self._log_value(gammalib.NORMAL, 'Append observation', '%s ("%s")' %
                            (obs_id, object_name))

            # Log the file names of the event and response files
            self._log_value(gammalib.EXPLICIT, ' Event file', eventfile)
            self._log_value(gammalib.EXPLICIT, ' Effective area file', aefffile)
            self._log_value(gammalib.EXPLICIT, ' Point spread function file',
                            psffile)
            self._log_value(gammalib.EXPLICIT, ' Energy dispersion file',
                            edispfile)
            self._log_value(gammalib.EXPLICIT, ' Background file', bkgfile)
            if self._use_bkg_scale:
                self._log_value(gammalib.EXPLICIT, ' Background scale',
                                bkg_scale)

        # Append models provided by "inmodel" to the model container
        self._append_inmodels()

        # Write observation summary
        self._write_summary()

        # Write warning in no observation could be used from runlist
        if self._xml_obslist.size() == 0:
            self._log_string(gammalib.NORMAL, 'WARNING: No observation was '
                             'appended from specified runlist')

        # Return
        return       
    
    def execute(self):
        """
        Execute the script
        """
        # Open logfile
        self.logFileOpen()

        # Read ahead output parameters
        self._read_ahead(True)

        # Run the script
        self.run()

        # Save residual map
        self.save()

        # Return
        return

    def save(self):
        """
        Save observation definition and model definition XML file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save XML files')

        # Get filenames
        outobs   = self['outobs'].filename()
        outmodel = self['outmodel'].filename()

        # Log observation definition file name
        self._log_value(gammalib.NORMAL, 'Obs. definition XML file',
                        outobs.url())

        # Save observation definition XML file
        self._xml.save(outobs.url())

        # Log model definition file name
        self._log_value(gammalib.NORMAL, 'Model definiton XML file',
                        outmodel.url())

        # Save model XML file
        self._models.save(outmodel.url())

        # Return
        return       

    def runlist(self, runlist):
        """
        Set runlist
        
        Parameters
        ----------
        runlist : list
            List of observation IDs

        Raises
        ------
        RuntimeError
            Input runlist is not a Python list
        """
        # If the runlist is convertable to a Python list then store the
        # runlist in the class
        try:
            self._runlist = list(runlist)

        # ... otherwise raise an exception
        except:
            raise RuntimeError('Argument is not a Python list')
        
        # Return
        return

    def ebounds(self):
        """
        Return runlist energy boundaries
        """
        # Return energy boundaries
        return self._ebounds
    
    def obs(self):
        """
        Return observations container
        """   
        # Clear observations
        self._obs.clear()
        
        # Read observations if XML is filled
        if self._xml.size():   
            self._obs.read(self._xml)
                
        # Assign models
        self._obs.models(self._models)
        
        # Return observations
        return self._obs


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csiactobs(sys.argv)

    # Execute application
    app.execute()
