#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csiactobs script
#
# Copyright (C) 2016 Juergen Knoedlseder
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
import cscripts
from testing import test


# =============================== #
# Test class for csiactobs script #
# =============================== #
class Test(test):
    """
    Test class for csiactobs script

    This test class makes unit tests for the csiactobs script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set data members
        self._datapath = self._datadir + '/iactdata'
        self._runlist  = self._datadir + '/iact_runlist.dat'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csiactobs')

        # Append tests
        self.append(self._test_cmd, 'Test csiactobs on command line')
        self.append(self._test_python, 'Test csiactobs from Python')

        # Return
        return

    # Test csiactobs on command line
    def _test_cmd(self):
        """
        Test csiactobs on the command line
        """
        # Set script name
        csiactobs = self._script('csiactobs')

        # Setup csiactobs command
        cmd = csiactobs+' datapath="'+self._datapath+'"'+ \
                        ' prodname="unit-test"'+ \
                        ' infile="'+self._runlist+'"'+ \
                        ' bkgpars=1'+\
                        ' outobs="csiactobs_obs_cmd1.xml"'+ \
                        ' outmodel="csiactobs_bgd_cmd1.xml"'+ \
                        ' logfile="csiactobs_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check observation definition XML file
        self._check_obsdef('csiactobs_obs_cmd1.xml', 6)

        # Check model definition XML file
        self._check_moddef('csiactobs_bgd_cmd1.xml', 6)

        # Setup csiactobs command
        cmd = csiactobs+' datapath="data_path_that_does_not_exist"'+ \
                        ' prodname="unit-test"'+ \
                        ' infile="'+self._runlist+'"'+ \
                        ' bkgpars=1'+\
                        ' outobs="csiactobs_obs_cmd2.xml"'+ \
                        ' outmodel="csiactobs_bgd_cmd2.xml"'+ \
                        ' logfile="csiactobs_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input datapath when executed from command line')
        
        # Setup csiactobs command
        cmd = csiactobs+' datapath="'+self._datapath+'"'+ \
                        ' prodname="unit-test-doesnt-exist"'+ \
                        ' infile="'+self._runlist+'"'+ \
                        ' bkgpars=1'+\
                        ' outobs="csiactobs_obs_cmd3.xml"'+ \
                        ' outmodel="csiactobs_bgd_cmd3.xml"'+ \
                        ' logfile="csiactobs_cmd3.log" chatter=1'       

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input prodname when executed from command line')

        # Return
        return

    # Test csiactobs from Python
    def _test_python(self):
        """
        Test csiactobs from Python
        """
        # Allocate empty csiactobs script
        iactobs = cscripts.csiactobs()

        # Check that empty csiactobs sciript has an empty observation container
        # and energy boundaries
        self.test_value(iactobs.obs().size(), 0,
             'Check that empty csiactobs has an empty observation container')
        self.test_value(iactobs.ebounds().size(), 0,
             'Check that empty csiactobs has empty energy bins')

        # Check that saving saves an empty model definition file
        iactobs['outobs']   = 'csiactobs_obs_py0.xml'
        iactobs['outmodel'] = 'csiactobs_bgd_py0.xml'
        iactobs['logfile']  = 'csiactobs_py0.log'
        iactobs.logFileOpen()
        iactobs.save()

        # Check empty observation definition XML file
        self._check_obsdef('csiactobs_obs_py0.xml', 0)

        # Check empty model definition XML file
        self._check_moddef('csiactobs_bgd_py0.xml', 0)

        # Check that clearing does not lead to an exception or segfault
        #iactobs.clear()

        # Set-up csiactobs
        iactobs = cscripts.csiactobs()
        iactobs['datapath'] = self._datapath
        iactobs['prodname'] = 'unit-test'
        iactobs['infile']   = self._runlist
        iactobs['bkgpars']  = 1
        iactobs['outobs']   = 'csiactobs_obs_py1.xml'
        iactobs['outmodel'] = 'csiactobs_bgd_py1.xml'
        iactobs['logfile']  = 'csiactobs_py1.log'
        iactobs['chatter']  = 2

        # Run csiactobs script and save run list
        iactobs.logFileOpen()   # Make sure we get a log file
        iactobs.run()
        iactobs.save()

        # Check observation definition XML file
        self._check_obsdef('csiactobs_obs_py1.xml', 6)
        
        # Check model definition XML file
        self._check_moddef('csiactobs_bgd_py1.xml', 6)
        
        # Create test runlist
        runlist = ['15000','15001']
           
        # Set-up csiactobs using a runlist with 2 background parameters
        iactobs = cscripts.csiactobs()
        iactobs['datapath'] = self._datapath
        iactobs['prodname'] = 'unit-test'
        iactobs['bkgpars']  = 2
        iactobs['outobs']   = 'csiactobs_obs_py2.xml'
        iactobs['outmodel'] = 'csiactobs_bgd_py2.xml'
        iactobs['logfile']  = 'csiactobs_py2.log'
        iactobs['chatter']  = 3
        iactobs.runlist(runlist)
   
        # Run csiactobs script and save run list
        iactobs.logFileOpen()   # Make sure we get a log file
        iactobs.run()
        iactobs.save()
           
        # Test return functions
        self.test_value(iactobs.obs().size(), 2,
                        'Check number of observations in container')
        self.test_value(iactobs.ebounds().size(), 0,
                        'Check number of energy boundaries')
   
        # Check observation definition XML file
        self._check_obsdef('csiactobs_obs_py2.xml',2)
           
        # Check model definition XML file
        self._check_moddef('csiactobs_bgd_py2.xml',2)
           
        # Set-up csiactobs with a large number of free parameters and "aeff"
        # background
        iactobs = cscripts.csiactobs()
        iactobs['datapath']      = self._datapath
        iactobs['prodname']      = 'unit-test'
        iactobs['infile']        = self._runlist
        iactobs['bkgpars']       = 8
        iactobs['bkg_mod_hiera'] = 'aeff'
        iactobs['outobs']        = 'csiactobs_obs_py3.xml'
        iactobs['outmodel']      = 'csiactobs_bgd_py3.xml'
        iactobs['logfile']       = 'csiactobs_py3.log'
        iactobs['chatter']       = 4
    
        # Execute csiactobs script
        iactobs.execute()
        
        # Check observation definition XML file
        self._check_obsdef('csiactobs_obs_py3.xml',6)
            
        # Check model definition XML file
        self._check_moddef('csiactobs_bgd_py3.xml',6)
            
        # Set-up csiactobs with a "gauss" background and "inmodel" parameter
        iactobs = cscripts.csiactobs()
        iactobs['datapath']      = self._datapath
        iactobs['inmodel']       = self._model
        iactobs['prodname']      = 'unit-test'
        iactobs['infile']        = self._runlist
        iactobs['bkgpars']       = 1
        iactobs['bkg_mod_hiera'] = 'gauss'
        iactobs['outobs']        = 'NONE'
        iactobs['outmodel']      = 'NONE'
        iactobs['logfile']       = 'csiactobs_py4.log'
        iactobs['chatter']       = 4
   
        # Run csiactobs script
        iactobs.logFileOpen()   # Make sure we get a log file
        iactobs.run()
        
        # Check number of observations
        self.test_value(iactobs.obs().size(), 6,
                        'Check number of observations in container')

        # Check number of models
        self.test_value(iactobs.obs().models().size(), 8,
                        'Check number of models in container')        
        
        # Set-up csiactobs with a "gauss" background and "inmodel" parameter
        iactobs = cscripts.csiactobs()
        iactobs['datapath']      = self._datapath
        iactobs['inmodel']       = self._model
        iactobs['prodname']      = 'unit-test'
        iactobs['infile']        = self._runlist
        iactobs['bkgpars']       = 1
        iactobs['bkg_mod_hiera'] = 'irf'
        iactobs['outobs']        = 'NONE'
        iactobs['outmodel']      = 'NONE'
        iactobs['logfile']       = 'csiactobs_py4.log'
        iactobs['chatter']       = 4
   
        # Run csiactobs script
        iactobs.logFileOpen()   # Make sure we get a log file
        iactobs.run()
        
        # Check number of observations
        self.test_value(iactobs.obs().size(), 5,
                        'Check number of observations in container')

        # Check number of models
        self.test_value(iactobs.obs().models().size(), 7,
                        'Check number of models in container')     

        # Return
        return

    # Check observation definition XML file
    def _check_obsdef(self, filename, obs_expected):
        """
        Check observation definition XML file
        """
        # Load observation definition XML file
        obs = gammalib.GObservations(filename)

        # Check number of observations
        self.test_value(obs.size(), obs_expected,
                        'Check for '+str(obs_expected)+' observations in XML file')
        
        # If there are observations in the XML file then check their content
        if obs_expected > 0:
            
            # Get response
            rsp = obs[0].response()
            
            # Test response
            self.test_value(obs[0].eventfile().file(), 'events_0.fits.gz',
                            'Check event file name')
            self.test_value(obs[0].eventfile().extname(), 'EVENTS',
                            'Check event extension name')
            self.test_value(rsp.aeff().filename().file(), 'irf_file.fits.gz',
                            'Check effective area file name')
            self.test_value(rsp.aeff().filename().extname(), 'EFFECTIVE AREA',
                            'Check effective area extension name')
            self.test_value(rsp.psf().filename().file(), 'irf_file.fits.gz',
                            'Check point spread function file name')
            self.test_value(rsp.psf().filename().extname(), 'POINT SPREAD FUNCTION',
                            'Check point spread function extension name')
            self.test_value(rsp.edisp().filename().file(), 'irf_file.fits.gz',
                            'Check energy dispersion file name')
            self.test_value(rsp.edisp().filename().extname(), 'ENERGY DISPERSION',
                            'Check energy dispersion extension name')
            self.test_value(rsp.background().filename().file(), 'irf_file.fits.gz',
                            'Check background file name')
            self.test_value(rsp.background().filename().extname(), 'BACKGROUND',
                            'Check background extension name')
        
        # Return
        return

    # Check model XML file
    def _check_moddef(self, filename, models_expected):
        """
        Check model definition XML file
        """
        # Load model definition XML file
        models = gammalib.GModels(filename)

        # Check number of models
        self.test_value(models.size(), models_expected,
                        'Check for '+str(models_expected)+' models in XML file')

        # Return
        return
