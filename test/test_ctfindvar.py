#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctfindvar tool.
#
# Copyright (C) 2018 Simon Bonnefoy
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
                        ' logfile="ctfindvar_cmd1.log" chatter=1'+\
                        ' prefix=ctfindvar_test_ '+\
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

        #setup ctfindvar
        ctfindvar['inobs']     = self._events
        ctfindvar['logfile']   = 'ctfindvar_py0.log'
        ctfindvar['chatter']   =  1
        ctfindvar['prefix']    = 'ctfindvar_test_python_'
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

        ## Run ctfindvar tool
        ctfindvar.logFileOpen()
        ctfindvar.run()

        # Check that saving does not nothing
        ctfindvar['logfile'] = 'ctfindvar_py1.log'
        ctfindvar.logFileOpen()
        ctfindvar.save()
 
        # Check result
        self._check_result_file('ctfindvar_test_python_srcsig.fits')

        # Copy ctfindvar tool
        cpy_tool = ctfindvar

        # Run copy of ctfindvar tool again
        cpy_tool['logfile'] = 'ctfindvar_py2.log'
        cpy_tool['chatter'] = 3
        cpy_tool['prefix']  = 'ctfindvar_test_copytool_' 
        cpy_tool.logFileOpen()
        cpy_tool.run()
        cpy_tool.save()

        # Check result
        self._check_result_file('ctfindvar_test_copytool_srcsig.fits')

        # Check that clearing does not lead to an exception or segfault
        ctfindvar.clear()

        # Return
        return

    # Check ctfindvar result
    def _check_result_file(self, filename):
        """
        Check content of tool

        Parameters
        ----------
        tool : `~ctools.ctfindvar`
            ctfindvar instance
        """
        print "checking if you amazing file is valid"
        # Read variability  file
        fits    = gammalib.GFits(filename)
        self.test_value(fits.size(), 3,
             'Check for 3 rows in srcsig file')
        # Return
        return
