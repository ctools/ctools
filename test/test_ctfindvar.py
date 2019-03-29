#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctfindvar tool.
#
# Copyright (C) 2018-2019 Simon Bonnefoy
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
import ctools
from testing import test


# ======================= #
# Test class for ctfindvar tool #
# ======================= #
class Test(test):
    """
    Test class for ctfindvar tool

    This test class makes unit tests for the ctfindvar tool by using it from
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
        self.name('ctfindvar')

        # Append tests
        self.append(self._test_cmd, 'Test ctfindvar on command line')
        self.append(self._test_python, 'Test ctfindvar from Python')

        # Return
        return

    # Test ctfindvar on command line
    def _test_cmd(self):
        """
        Test ctfindvar on the command line
        """
        # Set tool name
        ctfindvar = self._tool('ctfindvar')

        # Setup ctfindvar command
        cmd = ctfindvar+' inobs="'+self._events+'"'+\
                        ' outmap="ctfindvar_cmd1.fits"'+\
                        ' outmodel="ctfindvar_cmd1.xml"'+\
                        ' logfile="ctfindvar_cmd1.log" chatter=1'+\
                        ' tmin=NONE tmax=NONE tinterval=20'+\
                        ' emin=0.2 emax=1.1'+\
                        ' nxpix=200 nypix=200 binsz=0.05'+\
                        ' proj=CAR coordsys=CEL'+\
                        ' xref=83.63 yref=22.51'+\
                        ' xsrc=83.63 ysrc=22.51'+\
                        ' caldb=prod2 irf=South_0.5h'
        
        # Check if execution was successful
        self.test_value(self._execute(cmd), 0,
             'Check successful execution from command line')

        # Setup ctfindvar command
        cmd = ctfindvar+' inobs="event_file_that_does_not_exist.fits"'+ \
                        ' outmap="ctfindvar_cmd2.fits"'+ \
                        ' outmodel="ctfindvar_cmd2.xml"'+ \
                        ' logfile="ctfindvar_cmd2.log" chatter=1'+ \
                        ' tmin=NONE tmax=NONE tinterval=20'+ \
                        ' emin=0.2 emax=1.1'+ \
                        ' nxpix=200 nypix=200 binsz=0.05'+ \
                        ' proj=CAR coordsys=CEL'+ \
                        ' xref=83.63 yref=22.51'+ \
                        ' xsrc=83.63 ysrc=22.51'+ \
                        ' caldb=prod2 irf=South_0.5h'
        
        # Check if execution was successful
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check ctfindvar --help
        self._check_help(ctfindvar)

        # Return
        return

    # Test ctfindvar from Python
    def _test_python(self):
        """
        Test ctfindvar from Python
        """

        # Allocate ctfindvar
        ctfindvar = ctools.ctfindvar()

        # Setup ctfindvar
        ctfindvar['inobs']     = self._events
        ctfindvar['tmin']      = 'NONE'
        ctfindvar['tmax']      = 'NONE' 
        ctfindvar['tinterval'] = 20 
        ctfindvar['emin']      = 0.2 
        ctfindvar['emax']      = 1.1 
        ctfindvar['nxpix']     = 200 
        ctfindvar['nypix']     = 200 
        ctfindvar['binsz']     = 0.05 
        ctfindvar['coordsys']  = 'CEL'
        ctfindvar['proj']      = 'CAR'
        ctfindvar['xref']      = 83.63 
        ctfindvar['yref']      = 22.51
        ctfindvar['xsrc']      = 83.63
        ctfindvar['ysrc']      = 22.51
        ctfindvar['caldb']     = 'prod2'
        ctfindvar['irf']       = 'South_0.5h' 
        ctfindvar['outmap']    = 'ctfindvar_py0.fits'
        ctfindvar['outmodel']  = 'ctfindvar_py0.xml'
        ctfindvar['logfile']   = 'ctfindvar_py0.log'
        ctfindvar['chatter']   =  1

        # Run ctfindvar tool
        ctfindvar.logFileOpen()
        ctfindvar.run()

        # Check that saving does not nothing
        ctfindvar['outmap']   = 'ctfindvar_py1.fits'
        ctfindvar['outmodel'] = 'ctfindvar_py1.xml'
        ctfindvar['logfile']  = 'ctfindvar_py1.log'
        ctfindvar.logFileOpen()
        ctfindvar.save()
 
        # Check result
        self._check_result_file('ctfindvar_py1.fits')

        # Copy ctfindvar tool
        cpy_tool = ctfindvar.copy()

        # Run copy of ctfindvar tool again with lower threshold
        cpy_tool['threshold'] = 1.0
        cpy_tool['outmap']    = 'ctfindvar_py2.fits'
        cpy_tool['outmodel']  = 'ctfindvar_py2.xml'
        cpy_tool['logfile']   = 'ctfindvar_py2.log'
        cpy_tool['chatter']   = 3
        cpy_tool.logFileOpen()
        cpy_tool.execute()

        # Check result
        self._check_result_file('ctfindvar_py2.fits')

        # Check that clearing does not lead to an exception or segfault
        ctfindvar.clear()

        # Return
        return

    # Check ctfindvar result
    def _check_result_file(self, filename):
        """
        Check ctfindvar result

        Parameters
        ----------
        filename : str
            Output file
        """
        # Read variability file
        fits = gammalib.GFits(filename)
        self.test_value(fits.size(), 3, 'Check for 3 extensions in output file')

        # Return
        return
