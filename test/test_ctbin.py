#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctbin tool.
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


# ========================= #
# Test class for ctbin tool #
# ========================= #
class Test(test):
    """
    Test class for ctbin tool
    
    This test class makes unit tests for the ctbin tool by using it from
    the command line and from Python.
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
        self.name('ctbin')

        # Append tests
        self.append(self._test_cmd, 'Test ctbin on command line')
        self.append(self._test_python, 'Test ctbin from Python')

        # Return
        return

    # Test ctbin on command line
    def _test_cmd(self):
        """
        Test ctbin on the command line
        """
        # Set tool name
        ctbin = self._tool('ctbin')

        # Setup ctbin command
        cmd = ctbin+' inobs="'+self._events+'"'+ \
                    ' outcube="cntmap_cmd1.fits"'+\
                    ' emin=0.1 emax=100.0 enumbins=10 ebinalg="LOG"'+ \
                    ' nxpix=40 nypix=40 binsz=0.1 coordsys="CEL"'+ \
                    ' xref=83.63 yref=22.01 proj="CAR"'+ \
                    ' logfile="ctbin_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_value(self._execute(cmd), 0,
             'Check successful execution from command line')

        # Load counts cube and check content.
        evt = gammalib.GCTAEventCube('cntmap_cmd1.fits')
        self._check_cube(evt, 5542)

        # Setup ctbin command
        cmd = ctbin+' inobs="events_that_do_not_exist.fits"'+ \
                    ' outcube="cntmap_cmd2.fits"'+\
                    ' emin=0.1 emax=100.0 enumbins=10 ebinalg="LOG"'+ \
                    ' nxpix=40 nypix=40 binsz=0.1 coordsys="CEL"'+ \
                    ' xref=83.63 yref=22.01 proj="CAR"'+ \
                    ' logfile="ctbin_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctbin command with different ebounds than input
        cmd = ctbin+' inobs="'+self._events+'"'+ \
                    ' outcube="cntmap_cmd3.fits"'+\
                    ' emin=0.5 emax=20.0 enumbins=10 ebinalg="LOG"'+ \
                    ' nxpix=40 nypix=40 binsz=0.1 coordsys="CEL"'+ \
                    ' xref=83.63 yref=22.01 proj="CAR"'+ \
                    ' logfile="ctbin_cmd3.log" chatter=3'

        # Check if execution was successful
        self.test_value(self._execute(cmd), 0,
             'Check successful execution from command line')

        # Load counts cube and check content.
        evt = gammalib.GCTAEventCube('cntmap_cmd3.fits')
        self._check_cube(evt, 890)

        # Setup ctbin --help option
        cmd = ctbin+' --help'

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

    # Test ctbin from Python
    def _test_python(self):
        """
        Test ctbin from Python
        """
        # Allocate ctbin
        bin = ctools.ctbin()

        # Check that empty ctbin has an empty observation container and counts
        # cube
        self.test_value(bin.obs().size(), 0,
             'Check that empty ctbin has an empty observation container')
        self.test_value(bin.cube().size(), 0,
             'Check that empty ctbin has an empty counts cube')

        # Check that saving does not nothing
        bin['logfile'] = 'ctbin_py0.log'
        bin['outcube'] = 'ctbin_py0.fits'
        bin.logFileOpen()
        bin.save()
        self.test_assert(not os.path.isfile('ctbin_py0.fits'),
             'Check that no counts cube has been created')

        # Check that publishing does not lead to an exception or segfault
        bin.publish()

        # Check that clearing does not lead to an exception or segfault
        bin.clear()

        # Now set ctbin parameters
        bin['inobs']    = self._events
        bin['outcube']  = 'ctbin_py1.fits'
        bin['ebinalg']  = 'LOG'
        bin['emin']     = 0.1
        bin['emax']     = 100.0
        bin['enumbins'] = 10
        bin['nxpix']    = 40
        bin['nypix']    = 40
        bin['binsz']    = 0.1
        bin['coordsys'] = 'CEL'
        bin['proj']     = 'CAR'
        bin['xref']     = 83.63
        bin['yref']     = 22.01
        bin['logfile']  = 'ctbin_py1.log'
        bin['chatter']  = 2

        # Run ctbin tool
        bin.logFileOpen()
        bin.run()

        # Check content of observation container and counts cube
        self._check_observation(bin, 5542)
        self._check_cube(bin.cube(), 5542)

        # Copy ctbin tool
        cpy_bin = bin.copy()

        # Check content of observation container and counts cube
        self._check_observation(cpy_bin, 5542)
        self._check_cube(cpy_bin.cube(), 5542)

        # Run copy of ctbin tool again
        cpy_bin['logfile'] = 'ctbin_py2.log'
        cpy_bin['chatter'] = 3
        cpy_bin.logFileOpen()
        cpy_bin.run()

        # Check content of observation container and counts cube. Now an empty
        # event cube is expected since on input the observation is binned, and
        # any binned observation will be skipped.
        self._check_observation(cpy_bin, 0)
        self._check_cube(cpy_bin.cube(), 0)

        # Save counts cube
        bin.save()

        # Load counts cube and check content
        evt = gammalib.GCTAEventCube('ctbin_py1.fits')
        self._check_cube(evt, 5542)

        # Now clear ctbin tool
        bin.clear()

        # Check that cleared ctbin has an empty observation container and
        # counts cube
        self.test_value(bin.obs().size(), 0,
             'Check that empty ctbin has an empty observation container')
        self.test_value(bin.cube().size(), 0,
             'Check that empty ctbin has an empty counts cube')

        # Prepare observation container for stacking of events into a
        # single counts cube
        obs = self._obs_events()

        # Set-up ctbin using an observation container
        bin = ctools.ctbin(obs)
        bin['outcube']  = 'ctbin_py3.fits'
        bin['ebinalg']  = 'LOG'
        bin['emin']     = 0.1
        bin['emax']     = 100.0
        bin['enumbins'] = 10
        bin['nxpix']    = 40
        bin['nypix']    = 40
        bin['binsz']    = 0.1
        bin['coordsys'] = 'CEL'
        bin['proj']     = 'CAR'
        bin['xref']     = 83.63
        bin['yref']     = 22.01
        bin['publish']  = True
        bin['logfile']  = 'ctbin_py3.log'
        bin['chatter']  = 3

        # Execute ctbin tool
        bin.logFileOpen()
        bin.execute()

        # Check content of observation and cube (need multiplier=3 since
        # three identical observations have been appended)
        self._check_observation(bin, 5542, multiplier=3)
        self._check_cube(bin.cube(), 5542, multiplier=3)

        # Load counts cube and check content.
        evt = gammalib.GCTAEventCube('ctbin_py3.fits')
        self._check_cube(evt, 5542, multiplier=3)

        # And finally go for a fully Pythonic version with all parameters
        # being specified in a dictionary
        pars = {'inobs': self._events, 'ebinalg': 'LIN', 'emin': 0.1,
                'emax': 100.0, 'enumbins': 10, 'nxpix': 40, 'nypix': 40,
                'binsz': 0.1, 'coordsys': 'CEL', 'proj': 'CAR',
                'xref': 83.63, 'yref': 22.01, 'outcube': 'ctbin_py4.fits',
                'logfile': 'ctbin_py4.log', 'chatter': 2}
        bin = ctools.ctbin()
        bin.pardict(pars)
        bin.logFileOpen()
        bin.execute()

        # Load counts cube and check content
        evt = gammalib.GCTAEventCube('ctbin_py4.fits')
        self._check_cube(evt, 5542)

        # Return
        return

    # Check observation
    def _check_observation(self, ctbin, nevents, multiplier=1):
        """
        Check content of an observation
        
        Parameters
        ----------
        ctbin : `~ctools.ctbin`
            ctbin instance
        nevents : int
            Expected number of events
        multiplier : float, optional
            Ontime and livetime multiplier
        """
        # Test observation container
        obs = gammalib.GCTAObservation(ctbin.obs()[0])
        pnt = obs.pointing()
        self.test_value(ctbin.obs().size(), 1,
             'Check that there is one observation on the observation container')
        self.test_value(obs.instrument(), 'CTA',
             'Check that the one observation is a CTA observation')
        self.test_value(obs.ontime(), 1800.0*multiplier, 1.0e-6,
             'Check ontime of observation')
        self.test_value(obs.livetime(), 1710.0*multiplier, 1.0e-6,
             'Check livetime of observation')
        self.test_value(pnt.dir().ra_deg(), 83.63, 1.0e-6,
             'Check pointing Right Ascension of observation')
        self.test_value(pnt.dir().dec_deg(), 22.01, 1.0e-6,
             'Check pointing Declination of observation')

        # Test event cube
        self._check_cube(obs.events(), nevents, multiplier=multiplier)

        # Return
        return

    # Check event cube
    def _check_cube(self, cube, nevents, multiplier=1):
        """
        Check content of event cube
        
        Parameters
        ----------
        cube : `~gammalib.GCTAEventCube`
            Event cube
        nevents : int
            Expected number of events
        multiplier : float
            Event number multiplier
        """
        # Test event cube
        self.test_value(cube.size(), 16000,
             'Check bin size of counts cube')
        self.test_value(cube.number(), nevents*multiplier,
             'Check number of events in counts cube')

        # Return
        return
