#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctskymap tool.
#
# Copyright (C) 2014-2018 Juergen Knoedlseder
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
import ctools
from testing import test


# ============================ #
# Test class for ctskymap tool #
# ============================ #
class Test(test):
    """
    Test class for ctskymap tool
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set test data
        self._skymap_ring      = self._datadir + '/skymap_ring.fits'
        self._exclusion_region = self._datadir + '/exclusion.reg'
        self._exclusion_map    = self._datadir + '/exclusion.fits'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('ctskymap')

        # Append tests
        self.append(self._test_cmd, 'Test ctskymap on command line')
        self.append(self._test_python, 'Test ctskymap from Python')
        self.append(self._test_python_irf,
                    'Test ctskymap from Python for IRF background')
        self.append(self._test_python_ring,
                    'Test ctskymap from Python for RING background')
        
        # Return
        return

    # Test ctskymap on command line
    def _test_cmd(self):
        """
        Test ctskymap on the command line
        """
        # Set tool name
        ctskymap = self._tool('ctskymap')

        # Setup ctskymap command
        cmd = ctskymap+' inobs="'+self._events+'"'+\
                       ' outmap="ctskymap_cmd1.fits"'+ \
                       ' emin=0.1 emax=100.0 nxpix=20 nypix=20'+ \
                       ' binsz=0.2 coordsys="CEL" proj="CAR"'+ \
                       ' xref=83.63 yref=22.01 bkgsubtract="NONE"'+ \
                       ' logfile="ctskymap_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctskymap_cmd1.fits')

        # Setup ctskymap command
        cmd = ctskymap+' inobs="event_file_that_does_not_exist.fits"'+\
                       ' outmap="ctskymap_cmd2.fits"'+ \
                       ' emin=0.1 emax=100.0 nxpix=20 nypix=20'+ \
                       ' binsz=0.2 coordsys="CEL" proj="CAR"'+ \
                       ' xref=83.63 yref=22.01 bkgsubtract="NONE"'+ \
                       ' logfile="ctskymap_cmd2.log" debug=yes chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check ctskymap --help
        self._check_help(ctskymap)

        # Return
        return

    # Test ctskymap from Python
    def _test_python(self):
        """
        Test ctskymap from Python
        """
        # Allocate ctskymap
        skymap = ctools.ctskymap()

        # Check that empty ctskymap tool holds a map that has no pixels
        self._check_map(skymap.skymap(), nx=0, ny=0)

        # Check that saving does not nothing
        skymap['outmap']  = 'ctskymap_py0.fits'
        skymap['logfile'] = 'ctskymap_py0.log'
        skymap.logFileOpen()
        skymap.save()
        self.test_assert(not os.path.isfile('ctskymap_py0.fits'),
             'Check that no sky map has been created')

        # Check that clearing does not lead to an exception or segfault
        skymap.clear()

        # Now set ctskymap parameters
        skymap['inobs']       = self._events
        skymap['emin']        = 0.1
        skymap['emax']        = 100
        skymap['nxpix']       = 20
        skymap['nypix']       = 20
        skymap['binsz']       = 0.2
        skymap['coordsys']    = 'CEL'
        skymap['proj']        = 'CAR'
        skymap['xref']        = 83.63
        skymap['yref']        = 22.01
        skymap['bkgsubtract'] = 'NONE'
        skymap['outmap']      = 'ctskymap_py1.fits'
        skymap['logfile']     = 'ctskymap_py1.log'
        skymap['chatter']     = 2

        # Run ctskymap tool
        skymap.logFileOpen()   # Make sure we get a log file
        skymap.run()
        skymap.save()

        # Check result file
        self._check_result_file('ctskymap_py1.fits')

        # Copy ctskymap tool
        cpy_skymap = skymap.copy()

        # Check sky map of ctskymap copy
        self._check_map(cpy_skymap.skymap())

        # Execute copy of ctskymap tool again, now with a higher chatter
        # level than before
        cpy_skymap['usepnt']   = True
        cpy_skymap['emin']     = 0.2
        cpy_skymap['emax']     = 150.0
        cpy_skymap['coordsys'] = 'GAL'
        cpy_skymap['proj']     = 'CAR'
        cpy_skymap['outmap']   = 'ctskymap_py2.fits'
        cpy_skymap['logfile']  = 'ctskymap_py2.log'
        cpy_skymap['publish']  = True
        cpy_skymap['chatter']  = 3
        cpy_skymap.logFileOpen()  # Needed to get a new log file
        cpy_skymap.execute()

        # Check result file
        self._check_result_file('ctskymap_py2.fits')

        # Now clear copy of ctskymap tool
        cpy_skymap.clear()

        # Check that cleared ctskymap tool holds a map that has no pixels
        self._check_map(cpy_skymap.skymap(), nx=0, ny=0)

        # Get mixed observation container
        obs = self._obs_mixed()

        # Allocate ctskymap tool from observation container
        skymap = ctools.ctskymap(obs)
        skymap['emin']        = 0.1
        skymap['emax']        = 100
        skymap['nxpix']       = 20
        skymap['nypix']       = 20
        skymap['binsz']       = 0.2
        skymap['coordsys']    = 'GAL'
        skymap['proj']        = 'CAR'
        skymap['xref']        = 184.5575
        skymap['yref']        =  -5.7844
        skymap['bkgsubtract'] = 'NONE'
        skymap['outmap']      = 'ctskymap_py3.fits'
        skymap['logfile']     = 'ctskymap_py3.log'
        skymap['chatter']     = 4

        # Execute tool
        skymap.logFileOpen()  # Needed to get a new log file
        skymap.execute()

        # Check result file
        self._check_result_file('ctskymap_py3.fits')

        # Publish with name
        skymap.publish('My sky map')

        # Return
        return

    # Test ctskymap from Python for IRF background method
    def _test_python_irf(self):
        """
        Test ctskymap from Python for IRF background method
        """
        # Test IRF background subtraction
        skymap = ctools.ctskymap()
        skymap['inobs']       = self._events
        skymap['emin']        = 0.1
        skymap['emax']        = 100
        skymap['nxpix']       = 20
        skymap['nypix']       = 20
        skymap['binsz']       = 0.2
        skymap['coordsys']    = 'GAL'
        skymap['proj']        = 'CAR'
        skymap['xref']        = 184.5575
        skymap['yref']        =  -5.7844
        skymap['bkgsubtract'] = 'IRF'
        skymap['caldb']       = self._caldb
        skymap['irf']         = self._irf
        skymap['outmap']      = 'ctskymap_irf_py1.fits'
        skymap['logfile']     = 'ctskymap_irf_py1.log'
        skymap['chatter']     = 2

        # Execute tool
        skymap.logFileOpen()  # Needed to get a new log file
        skymap.execute()

        # Check result file
        self._check_result_file('ctskymap_irf_py1.fits')
        self._check_result_file('ctskymap_irf_py1.fits[BACKGROUND]')
        self._check_result_file('ctskymap_irf_py1.fits[SIGNIFICANCE]')

        # Return
        return

    # Test ctskymap from Python for RING background method
    def _test_python_ring(self):
        """
        Test ctskymap from Python for RING background method
        """
        # Test basic RING background subtraction
        skymap = ctools.ctskymap()
        skymap['inobs']       = self._events
        skymap['emin']        = 0.1
        skymap['emax']        = 100
        skymap['nxpix']       = 10
        skymap['nypix']       = 10
        skymap['binsz']       = 0.2
        skymap['coordsys']    = 'CEL'
        skymap['proj']        = 'CAR'
        skymap['xref']        = 83.63
        skymap['yref']        = 22.01
        skymap['bkgsubtract'] = 'RING'
        skymap['roiradius']   = 0.1
        skymap['inradius']    = 0.6
        skymap['outradius']   = 0.8
        skymap['iterations']  = 0
        skymap['inexclusion'] = 'NONE'
        skymap['caldb']       = self._caldb
        skymap['irf']         = self._irf
        skymap['outmap']      = 'ctskymap_ring_py1.fits'
        skymap['logfile']     = 'ctskymap_ring_py1.log'
        skymap['chatter']     = 2

        # Execute tool
        skymap.logFileOpen()  # Needed to get a new log file
        skymap.execute()

        # Check result file
        self._check_result_file('ctskymap_ring_py1.fits', 10, 10)
        self._check_result_file('ctskymap_ring_py1.fits[BACKGROUND]', 10, 10)
        self._check_result_file('ctskymap_ring_py1.fits[SIGNIFICANCE]', 10, 10)

        # Test iterative RING background subtraction
        skymap = ctools.ctskymap()
        skymap['inobs']       = self._events
        skymap['emin']        = 0.1
        skymap['emax']        = 100
        skymap['nxpix']       = 10
        skymap['nypix']       = 10
        skymap['binsz']       = 0.2
        skymap['coordsys']    = 'CEL'
        skymap['proj']        = 'CAR'
        skymap['xref']        = 83.63
        skymap['yref']        = 22.01
        skymap['bkgsubtract'] = 'RING'
        skymap['roiradius']   = 0.1
        skymap['inradius']    = 0.6
        skymap['outradius']   = 0.8
        skymap['iterations']  = 2
        skymap['threshold']   = 5.0
        skymap['inexclusion'] = 'NONE'
        skymap['caldb']       = self._caldb
        skymap['irf']         = self._irf
        skymap['outmap']      = 'ctskymap_ring_py2.fits'
        skymap['logfile']     = 'ctskymap_ring_py2.log'
        skymap['chatter']     = 2

        # Execute tool
        skymap.logFileOpen()  # Needed to get a new log file
        skymap.execute()

        # Check result file
        self._check_result_file('ctskymap_ring_py2.fits', 10, 10)
        self._check_result_file('ctskymap_ring_py2.fits[BACKGROUND]', 10, 10)
        self._check_result_file('ctskymap_ring_py2.fits[SIGNIFICANCE]', 10, 10)

        # Test iterative RING background subtraction, starting from a sky map
        skymap = ctools.ctskymap()
        skymap['inmap']       = self._skymap_ring
        skymap['bkgsubtract'] = 'RING'
        skymap['roiradius']   = 0.1
        skymap['inradius']    = 0.6
        skymap['outradius']   = 0.8
        skymap['iterations']  = 2
        skymap['threshold']   = 5.0
        skymap['inexclusion'] = 'NONE'
        skymap['outmap']      = 'ctskymap_ring_py3.fits'
        skymap['logfile']     = 'ctskymap_ring_py3.log'
        skymap['chatter']     = 2

        # Execute tool
        skymap.logFileOpen()  # Needed to get a new log file
        skymap.execute()

        # Check result file
        self._check_result_file('ctskymap_ring_py3.fits', 10, 10)
        self._check_result_file('ctskymap_ring_py3.fits[BACKGROUND]', 10, 10)
        self._check_result_file('ctskymap_ring_py3.fits[SIGNIFICANCE]', 10, 10)

        # Test iterative RING background subtraction with an exclusion region
        skymap = ctools.ctskymap()
        skymap['inmap']       = self._skymap_ring
        skymap['bkgsubtract'] = 'RING'
        skymap['roiradius']   = 0.1
        skymap['inradius']    = 0.6
        skymap['outradius']   = 0.8
        skymap['iterations']  = 2
        skymap['threshold']   = 5.0
        skymap['inexclusion'] = self._exclusion_region
        skymap['outmap']      = 'ctskymap_ring_py4.fits'
        skymap['logfile']     = 'ctskymap_ring_py4.log'
        skymap['chatter']     = 2

        # Execute tool
        skymap.logFileOpen()  # Needed to get a new log file
        skymap.execute()

        # Check result file
        self._check_result_file('ctskymap_ring_py4.fits', 10, 10)
        self._check_result_file('ctskymap_ring_py4.fits[BACKGROUND]', 10, 10)
        self._check_result_file('ctskymap_ring_py4.fits[SIGNIFICANCE]', 10, 10)

        # Test iterative RING background subtraction with an exclusion map
        skymap = ctools.ctskymap()
        skymap['inmap']       = self._skymap_ring
        skymap['bkgsubtract'] = 'RING'
        skymap['roiradius']   = 0.1
        skymap['inradius']    = 0.6
        skymap['outradius']   = 0.8
        skymap['iterations']  = 2
        skymap['threshold']   = 5.0
        skymap['inexclusion'] = self._exclusion_map
        skymap['outmap']      = 'ctskymap_ring_py5.fits'
        skymap['logfile']     = 'ctskymap_ring_py5.log'
        skymap['chatter']     = 2

        # Execute tool
        skymap.logFileOpen()  # Needed to get a new log file
        skymap.execute()

        # Check result file
        self._check_result_file('ctskymap_ring_py5.fits', 10, 10)
        self._check_result_file('ctskymap_ring_py5.fits[BACKGROUND]', 10, 10)
        self._check_result_file('ctskymap_ring_py5.fits[SIGNIFICANCE]', 10, 10)

        # Test RING background subtraction with direct significance computation
        skymap = ctools.ctskymap()
        skymap['inmap']       = self._skymap_ring
        skymap['bkgsubtract'] = 'RING'
        skymap['roiradius']   = 0.1
        skymap['inradius']    = 0.6
        skymap['outradius']   = 0.8
        skymap['iterations']  = 0
        skymap['threshold']   = 5.0
        skymap['inexclusion'] = 'NONE'
        skymap['usefft']      = False  # Direct computation, no FFT
        skymap['outmap']      = 'ctskymap_ring_py6.fits'
        skymap['logfile']     = 'ctskymap_ring_py6.log'
        skymap['chatter']     = 2

        # Execute tool
        skymap.logFileOpen()  # Needed to get a new log file
        skymap.execute()

        # Check result file
        self._check_result_file('ctskymap_ring_py6.fits', 10, 10)
        self._check_result_file('ctskymap_ring_py6.fits[BACKGROUND]', 10, 10)
        self._check_result_file('ctskymap_ring_py6.fits[SIGNIFICANCE]', 10, 10)

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nx=20, ny=20):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Sky map file name
        nx : int, optional
            Number of X pixels
        ny : int, optional
            Number of Y pixels
        """
        # Open result file
        skymap = gammalib.GSkyMap(filename)

        # Check sky map
        self._check_map(skymap, nx=nx, ny=ny)

        # Return
        return

    # Check sky map
    def _check_map(self, skymap, nx=20, ny=20):
        """
        Check sky map

        Parameters
        ----------
        skymap : `~gammalib.GSkyMap`
            Sky map
        nx : int, optional
            Number of X pixels
        ny : int, optional
            Number of Y pixels
        """
        # Determine number of maps
        if nx > 0 and ny > 0:
            nmaps = 1
        else:
            nmaps = 0
        
        # Check dimensions
        self.test_value(skymap.nmaps(), nmaps, 'Check number of maps')
        self.test_value(skymap.nx(), nx, 'Check for number of X pixels')
        self.test_value(skymap.ny(), ny, 'Check for number of Y pixels')

        # Return
        return
