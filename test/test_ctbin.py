#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctbin tool.
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

        # Set members
        self._inobs_two = self._datadir + '/obs_unbinned_two.xml'

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
                    ' outobs="cntmap_cmd1.fits"'+\
                    ' emin=1.0 emax=100.0 enumbins=10 ebinalg="LOG"'+ \
                    ' nxpix=40 nypix=40 binsz=0.1 coordsys="CEL"'+ \
                    ' xref=83.63 yref=22.01 proj="CAR"'+ \
                    ' logfile="ctbin_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_value(self._execute(cmd), 0,
             'Check successful execution from command line')

        # Load counts cube and check content.
        evt = gammalib.GCTAEventCube('cntmap_cmd1.fits')
        self._check_cube(evt, 245)

        # Setup ctbin command
        cmd = ctbin+' inobs="events_that_do_not_exist.fits"'+ \
                    ' outobs="cntmap_cmd2.fits"'+\
                    ' emin=1.0 emax=100.0 enumbins=10 ebinalg="LOG"'+ \
                    ' nxpix=40 nypix=40 binsz=0.1 coordsys="CEL"'+ \
                    ' xref=83.63 yref=22.01 proj="CAR"'+ \
                    ' logfile="ctbin_cmd2.log" debug=yes chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Setup ctbin command with different ebounds than input
        cmd = ctbin+' inobs="'+self._events+'"'+ \
                    ' outobs="cntmap_cmd3.fits"'+\
                    ' emin=2.0 emax=20.0 enumbins=10 ebinalg="LOG"'+ \
                    ' nxpix=40 nypix=40 binsz=0.1 coordsys="CEL"'+ \
                    ' xref=83.63 yref=22.01 proj="CAR"'+ \
                    ' logfile="ctbin_cmd3.log" chatter=3'

        # Check if execution was successful
        self.test_value(self._execute(cmd), 0,
             'Check successful execution from command line')

        # Load counts cube and check content.
        evt = gammalib.GCTAEventCube('cntmap_cmd3.fits')
        self._check_cube(evt, 115)

        # Check joint (unstacked) binning

        # Setup ctbin command
        cmd = ctbin+' inobs="'+self._inobs_two+'"'+ \
                    ' outobs="obs_unbinned_two_binned.xml"'+\
                    ' emin=1.0 emax=100.0 enumbins=10 ebinalg="LOG"'+ \
                    ' nxpix=40 nypix=40 binsz=0.1 coordsys="CEL"'+ \
                    ' xref=83.63 yref=22.01 proj="CAR"'+ \
                    ' logfile="ctbin_cmd4.log" chatter=1 stack=no'+ \
                    ' prefix="cntmap_cmd4_"'

        # Check if execution was successful
        self.test_value(self._execute(cmd), 0,
             'Check successful execution from command line with stack=no')

        # Load counts cubes and check content.
        evt = gammalib.GCTAEventCube('cntmap_cmd4_cta_00001.fits')
        self._check_cube(evt, 245)

        # Load counts cubes and check content.
        evt = gammalib.GCTAEventCube('cntmap_cmd4_cta_00002.fits')
        self._check_cube(evt, 245)

        # Check observation definition XML output file
        obs = gammalib.GObservations('obs_unbinned_two_binned.xml')
        self.test_value(obs.size(), 2, 'Check for 2 observations in XML file')

        # Check for content of observations file
        self.test_value(obs[0].eventfile().file(), 'cntmap_cmd4_cta_00001.fits',
                        'Check first counts cube file name')
        self.test_value(obs[1].eventfile().file(), 'cntmap_cmd4_cta_00002.fits',
                        'Check second counts cube file name')


        # Check ctbin --help
        self._check_help(ctbin)

        # Return
        return

    # Test ctbin from Python
    def _test_python(self):
        """
        Test ctbin from Python
        """
        # Allocate ctbin
        binning = ctools.ctbin()

        # Check that empty ctbin has an empty observation container and counts
        # cube
        self.test_value(binning.obs().size(), 0,
             'Check that empty ctbin has an empty observation container')
        self.test_value(binning.cubes(), 0,
             'Check that empty ctbin has an empty counts cube')

        # Check that saving does not nothing
        binning['logfile'] = 'ctbin_py0.log'
        binning['outobs']  = 'ctbin_py0.fits'
        binning.logFileOpen()
        binning.save()
        self.test_assert(not os.path.isfile('ctbin_py0.fits'),
             'Check that no counts cube has been created')

        # Check that publishing does not lead to an exception or segfault
        binning.publish()

        # Check that clearing does not lead to an exception or segfault
        binning.clear()

        # Now set ctbin parameters
        binning['inobs']    = self._events
        binning['outobs']   = 'ctbin_py1.fits'
        binning['ebinalg']  = 'LOG'
        binning['emin']     = 1.0
        binning['emax']     = 100.0
        binning['enumbins'] = 10
        binning['nxpix']    = 40
        binning['nypix']    = 40
        binning['binsz']    = 0.1
        binning['coordsys'] = 'CEL'
        binning['proj']     = 'CAR'
        binning['xref']     = 83.63
        binning['yref']     = 22.01
        binning['logfile']  = 'ctbin_py1.log'
        binning['chatter']  = 2

        # Run ctbin tool
        binning.logFileOpen()
        binning.run()

        # Check content of observation container and counts cube
        self._check_observation(binning, 245)
        self._check_cube(binning.cube(), 245)

        # Copy ctbin tool
        cpy_bin = binning.copy()

        # Check content of observation container and counts cube
        self._check_observation(cpy_bin, 245)
        self._check_cube(cpy_bin.cube(), 245)

        # Run copy of ctbin tool again
        cpy_bin['logfile'] = 'ctbin_py2.log'
        cpy_bin['chatter'] = 3
        cpy_bin.logFileOpen()
        cpy_bin.run()

        # Check content of observation container and number of counts cubes.
        # There should be a single binned observation in the observation
        # container, which is the one that was produced in the run before.
        # Since ctbin finds no unbinned observation in the container, the
        # number of cubes should be zero.
        self._check_observation(cpy_bin, 245)
        self.test_value(cpy_bin.cubes(), 0, 'Check that there are no counts cubes')

        # Save counts cube
        binning.save()

        # Load counts cube and check content
        evt = gammalib.GCTAEventCube('ctbin_py1.fits')
        self._check_cube(evt, 245)

        # Now clear ctbin tool
        binning.clear()

        # Check that cleared ctbin has an empty observation container and
        # counts cube
        self.test_value(binning.obs().size(), 0,
             'Check that empty ctbin has an empty observation container')
        self.test_value(binning.cubes(), 0,
             'Check that empty ctbin has an empty counts cube')

        # Prepare observation container for stacking of events into a
        # single counts cube
        obs = self._obs_events()

        # Set-up ctbin using an observation container
        binning = ctools.ctbin(obs)
        binning['outobs']   = 'ctbin_py3.fits'
        binning['ebinalg']  = 'LOG'
        binning['emin']     = 1.0
        binning['emax']     = 100.0
        binning['enumbins'] = 10
        binning['nxpix']    = 40
        binning['nypix']    = 40
        binning['binsz']    = 0.1
        binning['coordsys'] = 'CEL'
        binning['proj']     = 'CAR'
        binning['xref']     = 83.63
        binning['yref']     = 22.01
        binning['publish']  = True
        binning['logfile']  = 'ctbin_py3.log'
        binning['chatter']  = 3

        # Execute ctbin tool
        binning.logFileOpen()
        binning.execute()

        # Check content of observation and cube (need multiplier=3 since
        # three identical observations have been appended)
        self._check_observation(binning, 245, multiplier=3)
        self._check_cube(binning.cube(), 245, multiplier=3)

        # Load counts cube and check content.
        evt = gammalib.GCTAEventCube('ctbin_py3.fits')
        self._check_cube(evt, 245, multiplier=3)

        # Now go for a fully Pythonic version with all parameters being
        # specified in a dictionary
        pars = {'inobs': self._events, 'ebinalg': 'LIN', 'emin': 1.0,
                'emax': 100.0, 'enumbins': 10, 'nxpix': 40, 'nypix': 40,
                'binsz': 0.1, 'coordsys': 'CEL', 'proj': 'CAR',
                'xref': 83.63, 'yref': 22.01, 'outobs': 'ctbin_py4.fits',
                'logfile': 'ctbin_py4.log', 'chatter': 2}
        binning = ctools.ctbin()
        binning.pardict(pars)
        binning.logFileOpen()
        binning.execute()

        # Load counts cube and check content
        evt = gammalib.GCTAEventCube('ctbin_py4.fits')
        self._check_cube(evt, 245)

        # Test unstacked version
        binning = ctools.ctbin()
        binning['inobs']    = self._inobs_two
        binning['stack']    = False
        binning['prefix']   = 'cntcube_py5_'
        binning['ebinalg']  = 'LOG'
        binning['emin']     = 1.0
        binning['emax']     = 100.0
        binning['enumbins'] = 10
        binning['nxpix']    = 40
        binning['nypix']    = 40
        binning['binsz']    = 0.1
        binning['coordsys'] = 'CEL'
        binning['proj']     = 'CAR'
        binning['xref']     = 83.63
        binning['yref']     = 22.01
        binning['outobs']   = 'ctbin_py5.xml'
        binning['logfile']  = 'ctbin_py5.log'
        binning['chatter']  = 2
        binning.logFileOpen()
        binning.run()

        # Check individual cubes
        self._check_cube(binning.cube(0), 245)
        self._check_cube(binning.cube(1), 245)

        # Save counts cube
        binning.save()

        # Load observations
        obs = gammalib.GObservations('ctbin_py5.xml')
        self.test_value(obs.size(), 2, 'Check number of output observations')

        # Check counts cubes
        evt = obs[0].events()
        self._check_cube(evt, 245)
        evt = obs[1].events()
        self._check_cube(evt, 245)

        # Load counts cubes and check content
        evt = gammalib.GCTAEventCube('cntcube_py5_cta_00001.fits')
        self._check_cube(evt, 245)
        evt = gammalib.GCTAEventCube('cntcube_py5_cta_00002.fits')
        self._check_cube(evt, 245)

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
        self.test_value(obs.ontime(), 300.0*multiplier, 1.0e-6,
             'Check ontime of observation')
        self.test_value(obs.livetime(), 294.0*multiplier, 1.0e-6,
             'Check livetime of observation')
        self.test_value(pnt.dir().ra_deg(), 83.63, 1.0e-6,
             'Check pointing Right Ascension of observation')
        self.test_value(pnt.dir().dec_deg(), 22.51, 1.0e-6,
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
