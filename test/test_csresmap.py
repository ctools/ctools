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
import gammalib
import cscripts
from testing import test


# ============================== #
# Test class for csresmap script #
# ============================== #
class Test(test):
    """
    Test class for csresmap script

    This test class makes unit tests for the csresmap script by using it
    from the command line and from Python.
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
        # Set script name
        csresmap = self._script('csresmap')

        # Setup csresmap command
        cmd = csresmap+' inobs="'+self._events+'"'+ \
                       ' outmap="csresmap_cmd1.fits"'+ \
                       ' inmodel="'+self._model+'"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                       ' emin=0.1 emax=100.0 enumbins=20'+ \
                       ' nxpix=50 nypix=50 binsz=0.02'+ \
                       ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                       ' algorithm="SUBDIVSQRT"'+ \
                       ' logfile="csresmap_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('csresmap_cmd1.fits')

        # Setup csresmap command
        cmd = csresmap+' inobs="event_file_that_does_not_exist.fits"'+ \
                       ' outmap="csresmap_cmd1.fits"'+ \
                       ' inmodel="'+self._model+'"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                       ' emin=0.1 emax=100.0 enumbins=20'+ \
                       ' nxpix=50 nypix=50 binsz=0.02'+ \
                       ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                       ' algorithm="SUBDIVSQRT"'+ \
                       ' logfile="csresmap_cmd2.log" chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csresmap from Python
    def _test_python(self):
        """
        Test csresmap from Python
        """
        # Set-up csresmap for event list
        resmap = cscripts.csresmap()
        resmap['inobs']     = self._events
        resmap['outmap']    = 'csresmap_py1.fits'
        resmap['inmodel']   = self._model
        resmap['caldb']     = self._caldb
        resmap['irf']       = self._irf
        resmap['emin']      = 0.1
        resmap['emax']      = 100.0
        resmap['enumbins']  = 20
        resmap['nxpix']     = 50
        resmap['nypix']     = 50
        resmap['binsz']     = 0.02
        resmap['coordsys']  = 'CEL'
        resmap['proj']      = 'CAR'
        resmap['xref']      = 83.63
        resmap['yref']      = 22.01
        resmap['algorithm'] = 'SUB'
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
        resmap['inobs']     = self._cntcube
        resmap['modcube']   = 'NONE'
        resmap['expcube']   = 'NONE'
        resmap['psfcube']   = 'NONE'
        resmap['edispcube'] = 'NONE'
        resmap['bkgcube']   = 'NONE'
        resmap['caldb']     = self._caldb
        resmap['irf']       = self._irf
        resmap['inmodel']   = self._model
        resmap['outmap']    = 'csresmap_py2.fits'
        resmap['algorithm'] = 'SUBDIV'
        resmap['logfile']   = 'csresmap_py2.log'
        resmap['chatter']   = 3
        resmap['publish']   = True

        # Run csresmap script
        resmap.execute()

        # Check pull distribution file
        self._check_result_file('csresmap_py2.fits', nx=200, ny=200)

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nx=50, ny=50):
        """
        Check result file
        """
        # Load residual map
        residual = gammalib.GSkyMap(filename)

        # Check residual map
        self.test_value(residual.nmaps(), 1, 'One map')
        self.test_value(residual.nx(), nx, '%d pixels in X' % nx)
        self.test_value(residual.ny(), ny, '%d pixels in Y' % ny)

        # Return
        return
