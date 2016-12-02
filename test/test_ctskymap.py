#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctskymap tool.
#
# Copyright (C) 2014-2016 Juergen Knoedlseder
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
                       ' emin=0.1 emax=100.0 nxpix=200 nypix=200'+ \
                       ' binsz=0.02 coordsys="CEL" proj="CAR"'+ \
                       ' xref=83.63 yref=22.01 bkgsubtract="NONE"'+ \
                       ' logfile="ctskymap_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctskymap_cmd1.fits')

        # Setup ctskymap command
        cmd = ctskymap+' inobs="event_file_that_does_not_exist.fits"'+\
                       ' outmap="ctskymap_cmd2.fits"'+ \
                       ' emin=0.1 emax=100.0 nxpix=200 nypix=200'+ \
                       ' binsz=0.02 coordsys="CEL" proj="CAR"'+ \
                       ' xref=83.63 yref=22.01 bkgsubtract="NONE"'+ \
                       ' logfile="ctskymap_cmd2.log" chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctskymap --help option
        cmd = ctskymap+' --help'

        # Check if execution was successful in case that the CTOOLS
        # environment variable was set or failed otherwise
        if 'CTOOLS' in os.environ:
            self.test_value(self._execute(cmd), 0,
                 'Check successful execution with --help option')
        else:
            self.test_assert(self._execute(cmd) != 0,
                 'Check execution failure with --help option')

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
        skymap['nxpix']       = 200
        skymap['nypix']       = 200
        skymap['binsz']       = 0.02
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
        skymap['nxpix']       = 200
        skymap['nypix']       = 200
        skymap['binsz']       = 0.02
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

        # Allocate ctskymap tool from observation container
        skymap = ctools.ctskymap(obs)
        skymap['emin']        = 0.1
        skymap['emax']        = 100
        skymap['nxpix']       = 200
        skymap['nypix']       = 200
        skymap['binsz']       = 0.02
        skymap['coordsys']    = 'GAL'
        skymap['proj']        = 'CAR'
        skymap['xref']        = 184.5575
        skymap['yref']        =  -5.7844
        skymap['bkgsubtract'] =  'IRF'
        skymap['caldb']       = self._caldb
        skymap['irf']         = self._irf
        skymap['outmap']      = 'ctskymap_py4.fits'
        skymap['logfile']     = 'ctskymap_py4.log'
        skymap['chatter']     = 4

        # Execute tool
        skymap.logFileOpen()  # Needed to get a new log file
        skymap.execute()

        # Check result file
        self._check_result_file('ctskymap_py3.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nx=200, ny=200):
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
    def _check_map(self, skymap, nx=200, ny=200):
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
