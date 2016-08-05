#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csiactcopy script.
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
import os
import gammalib
import cscripts
from testing import test


# ================================ #
# Test class for csiactcopy script #
# ================================ #
class Test(test):
    """
    Test class for csiactcopy script

    This test class makes unit tests for the csiactcopy script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        test.__init__(self)

        # Set data members
        self._remote   = self._datadir + '/iactdata/master.json'
        self._runlist1 = self._datadir + '/iact_runlist_1.dat'
        self._runlist2 = self._datadir + '/iact_runlist_2.dat'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csiactcopy')

        # Append tests
        self.append(self._test_cmd, 'Test csiactcopy on command line')
        self.append(self._test_python, 'Test csiactcopy from Python')

        # Return
        return

    # Test csiactcopy on command line
    def _test_cmd(self):
        """
        Test csiactcopy on the command line
        """
        # Set script name
        csiactcopy = self._script('csiactcopy')

        # Setup csiactcopy command
        cmd = csiactcopy+' remote_master="'+self._remote+'"'+ \
                         ' prodname="unit-test"'+ \
                         ' outpath="csiactcopy_cmd1"'+ \
                         ' logfile="csiactcopy_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check copy
        self._check_copy('csiactcopy_cmd1')

        # Setup csiactcopy command
        cmd = csiactcopy+' remote_master="master_that_does_not_exist"'+ \
                         ' prodname="unit-test"'+ \
                         ' outpath="csiactcopy_cmd2"'+ \
                         ' logfile="csiactcopy_cmd2.log"'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csiactcopy from Python
    def _test_python(self):
        """
        Test csiactcopy from Python
        """
        # Set-up csiactcopy
        iactcopy = cscripts.csiactcopy()
        iactcopy['remote_master'] = self._remote
        iactcopy['prodname']      = 'unit-test'
        iactcopy['outpath']       = 'csiactcopy_py1'
        iactcopy['logfile']       = 'csiactcopy_py1.log'
        iactcopy['chatter']       = 2

        # Run csiactcopy script and save run list
        iactcopy.logFileOpen()   # Make sure we get a log file
        iactcopy.run()
        iactcopy.save()
        
        # Check IACT data store
        self._check_copy('csiactcopy_py1')
        self._check_n_obs('csiactcopy_py1', 9)
        
        # Set-up csiactcopy with clobber=True, should lead to skipping
        # of files
        iactcopy = cscripts.csiactcopy()
        iactcopy['remote_master'] = self._remote
        iactcopy['prodname']      = 'unit-test'
        iactcopy['outpath']       = 'csiactcopy_py1'
        iactcopy['logfile']       = 'csiactcopy_py2.log'
        iactcopy['chatter']       = 3
        iactcopy['clobber']       = True

        # Execute csiactcopy script
        iactcopy.logFileOpen()   # Make sure we get a log file
        iactcopy.execute()
        
        # Check copy
        self._check_copy('csiactcopy_py1')
        self._check_n_obs('csiactcopy_py1', 9)
        
        # Set-up csiactcopy with clobber=False, should also lead to skipping
        # of files, but for the reason that it's not allowed to overwrite
        # the file
        iactcopy = cscripts.csiactcopy()
        iactcopy['remote_master'] = self._remote
        iactcopy['prodname']      = 'unit-test'
        iactcopy['outpath']       = 'csiactcopy_py1'
        iactcopy['logfile']       = 'csiactcopy_py3.log'
        iactcopy['chatter']       = 4
        iactcopy['clobber']       = False

        # Run csiactcopy script and save run list
        iactcopy.logFileOpen()   # Make sure we get a log file
        iactcopy.run()
        iactcopy.save()
        
        # Check copy
        self._check_copy('csiactcopy_py1')
        self._check_n_obs('csiactcopy_py1', 9)
        
        # Remove the copied IACT data store before continuing
        if os.path.isdir('csiactcopy_py2'):
            os.system('rm -r csiactcopy_py2')
        
        # Set-up csiactcopy using a run list with 5 runs
        iactcopy = cscripts.csiactcopy()
        iactcopy['runlist']       = self._runlist1
        iactcopy['remote_master'] = self._remote
        iactcopy['prodname']      = 'unit-test'
        iactcopy['outpath']       = 'csiactcopy_py2'
        iactcopy['logfile']       = 'csiactcopy_py4.log'
        iactcopy['chatter']       = 4

        # Run csiactcopy script and save run list
        iactcopy.logFileOpen()   # Make sure we get a log file
        iactcopy.run()
        iactcopy.save()

        # Check copy
        self._check_copy('csiactcopy_py2')
        self._check_n_obs('csiactcopy_py2', 5)

        # Set-up csiactcopy using another run list with 4 runs. Should have
        # 9 runs in total in the data store afterwards
        iactcopy = cscripts.csiactcopy()
        iactcopy['runlist']       = self._runlist2
        iactcopy['remote_master'] = self._remote
        iactcopy['prodname']      = 'unit-test'
        iactcopy['outpath']       = 'csiactcopy_py2'
        iactcopy['logfile']       = 'csiactcopy_py5.log'
        iactcopy['chatter']       = 4
        
        # Run csiactcopy script and save run list
        iactcopy.logFileOpen()   # Make sure we get a log file
        iactcopy.run()
        iactcopy.save()

        # Check copy
        self._check_copy('csiactcopy_py2')
        self._check_n_obs('csiactcopy_py2', 9)
        
        # Set-up csiactcopy using the same run list with 4 runs again, should
        # lead to skipping of files, and we should still have 9 runs at the
        # end in the datastore.
        iactcopy = cscripts.csiactcopy()
        iactcopy['runlist']       = self._runlist2
        iactcopy['remote_master'] = self._remote
        iactcopy['prodname']      = 'unit-test'
        iactcopy['outpath']       = 'csiactcopy_py2'
        iactcopy['logfile']       = 'csiactcopy_py6.log'
        iactcopy['chatter']       = 4
        iactcopy['clobber']       = True
        
        # Run csiactcopy script and save run list
        iactcopy.logFileOpen()   # Make sure we get a log file
        iactcopy.run()
        iactcopy.save()

        # Check copy
        self._check_copy('csiactcopy_py2')
        self._check_n_obs('csiactcopy_py2',9)

        # Return
        return

    # Check copy of IACT data store
    def _check_copy(self, pathname):
        """
        Check copy of IACT data store

        Parameters
        ----------
        pathname : str
            Path to copied IACT data store
        """
        # Set file names
        hdu_index_name = gammalib.GFilename(pathname+'/hdu-index.fits')
        obs_index_name = gammalib.GFilename(pathname+'/obs-index.fits')
        master_name    = gammalib.GFilename(pathname+'/master.json')

        # Check for existence of files
        self.test_assert(hdu_index_name.exists(),
                         'Check if file "hdu-index.fits" exists')
        self.test_assert(obs_index_name.exists(),
                         'Check if file "obs-index.fits" exists')
        self.test_assert(master_name.exists(),
                         'Check if file "master.json" exists')

        # Check if index files are FITS file
        self.test_assert(hdu_index_name.is_fits(),
                         'Check if file "hdu-index.fits" is FITS file')
        self.test_assert(obs_index_name.is_fits(),
                         'Check if file "obs-index.fits" is FITS file')
        
        # Return
        return

    # Check copy
    def _check_n_obs(self, pathname, n_expected):
        """
        Check number of available observations

        Parameters
        ----------
        pathname : str
            Path to copied IACT data store
        n_expected : int
            Expected number of observations in IACT data store
        """
        # Set file name
        obs_index_name = gammalib.GFilename(pathname+'/obs-index.fits[OBS_INDEX]')
        
        # Open index file
        fits = gammalib.GFits(obs_index_name)
        
        # Get number of observations
        n_obs = fits[obs_index_name.extname()].nrows()
        
        # Close FITS file
        fits.close()

        # Check for existence of observations
        self.test_value(n_obs, n_expected, 'Check for number of observations')
                
        # Return
        return