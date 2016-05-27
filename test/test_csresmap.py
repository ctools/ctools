#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csresmap script.
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


# ============================== #
# Test class for csresmap script #
# ============================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for csresmap script.

    This test class makes unit tests for the csresmap script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('csresmap')

        # Append tests
        self.append(self._test_cmd, 'Test csresmap on command line')
        self.append(self._test_python, 'Test csresmap from Python')

        # Return
        return

    # Test csresmap on command line
    def _test_cmd(self):
        """
        Test csresmap on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile('README'):
            csresmap = '../cscripts/csresmap.py'
        else:
            csresmap = 'csresmap'

        # Setup csresmap command
        cmd = csresmap+' inobs="data/crab_events.fits"'+ \
                       ' outmap="csresmap_cmd1.fits"'+ \
                       ' inmodel="data/crab.xml"'+ \
                       ' caldb="prod2" irf="South_0.5h"' + \
                       ' nxpix=50 nypix=50 binsz=0.02'+ \
                       ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                       ' algorithm="SUBDIV"'+ \
                       ' logfile="csresmap_cmd1.log" chatter=1'

        # Execute csresmap, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         'Successful csresmap execution on command line')

        # Check result file
        self._check_result_file('csresmap_cmd1.fits')

        # Setup csresmap command
        cmd = csresmap+' inobs="event_file_that_does_not_exist.fits"'+ \
                       ' outmap="csresmap_cmd1.fits"'+ \
                       ' inmodel="data/crab.xml"'+ \
                       ' caldb="prod2" irf="South_0.5h"' + \
                       ' nxpix=50 nypix=50 binsz=0.02'+ \
                       ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                       ' algorithm="SUBDIV"'+ \
                       ' logfile="csresmap_cmd2.log" chatter=2'

        # Execute csresmap, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         'Failure of csresmap execution on command line')

        # Return
        return

    # Test csresmap from Python
    def _test_python(self):
        """
        Test csresmap from Python.
        """
        # Set-up csresmap for event list
        resmap = cscripts.csresmap()
        resmap['inobs']     = 'data/crab_events.fits'
        resmap['outmap']    = 'csresmap_py1.fits'
        resmap['inmodel']   = 'data/crab.xml'
        resmap['caldb']     = 'prod2'
        resmap['irf']       = 'South_0.5h'
        resmap['nxpix']     = 50
        resmap['nypix']     = 50
        resmap['binsz']     = 0.02
        resmap['coordsys']  = 'CEL'
        resmap['proj']      = 'CAR'
        resmap['xref']      = 83.63
        resmap['yref']      = 22.01
        resmap['algorithm'] = 'SUBDIV'
        resmap['logfile']   = 'csresmap_py1.log'
        resmap['chatter']   = 2

        # Run csresmap script
        resmap.logFileOpen()   # Make sure we get a log file
        resmap.run()
        resmap.save()

        # Check pull distribution file
        self._check_result_file('csresmap_py1.fits')

        # Set-up csresmap for counts cube
        resmap = cscripts.csresmap()
        resmap['inobs']     = 'data/crab_cntmap.fits'
        resmap['modcube']   = 'NONE'
        resmap['expcube']   = 'NONE'
        resmap['psfcube']   = 'NONE'
        resmap['edispcube'] = 'NONE'
        resmap['bkgcube']   = 'NONE'
        resmap['caldb']     = 'prod2'
        resmap['irf']       = 'South_0.5h'
        resmap['inmodel']   = 'data/crab.xml'
        resmap['outmap']    = 'csresmap_py2.fits'
        resmap['algorithm'] = 'SUBDIV'
        resmap['logfile']   = 'csresmap_py2.log'
        resmap['chatter']   = 3

        # Run csresmap script
        resmap.execute()

        # Check pull distribution file
        self._check_result_file('csresmap_py2.fits', nx=200, ny=200)

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nx=50, ny=50):
        """
        Check result file.
        """
        # Load residual map
        residual = gammalib.GSkyMap(filename)

        # Check residual map
        self.test_value(residual.nmaps(), 1, 'One map')
        self.test_value(residual.nx(), nx, '%d pixels in X' % nx)
        self.test_value(residual.ny(), ny, '%d pixels in Y' % ny)

        # Return
        return
