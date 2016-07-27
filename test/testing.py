# ==========================================================================
# Python base class for unit tests
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


# ================================ #
# Python base class for unit tests #
# ================================ #
class test(gammalib.GPythonTestSuite):
    """
    Python base class for unit tests
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Set test data directory
        self._datadir = os.environ['TEST_DATA']

        # Set some standard test data (additional test data will be set in
        # the derived classes as needed)
        self._events  = self._datadir + '/crab_events.fits'
        self._cntcube = self._datadir + '/crab_cntmap.fits'
        self._model   = self._datadir + '/crab.xml'
        self._caldb   = 'prod2'
        self._irf     = 'South_0.5h'

        # Return
        return

    # Set script name
    def _script(self, name):
        """
        Set script name
        """
        # Set script name dependent on whether the software is installed or
        # not. If the software is not installed the 'TEST_BUILDDIR' environment
        # variable will exist
        if 'TEST_BUILDDIR' in os.environ:
            script = os.environ['TEST_BUILDDIR'] + '/cscripts/'+name+'.py'
        else:
            script = name

        # Return script name
        return script

    # Set tool name
    def _tool(self, name):
        """
        Set tool name
        """
        # Set tool name dependent on whether the software is installed or
        # not. If the software is not installed the 'TEST_BUILDDIR' environment
        # variable will exist
        if 'TEST_BUILDDIR' in os.environ:
            tool = os.environ['TEST_BUILDDIR'] + '/src/'+name+'/'+name
        else:
            tool = name

        # Return tool name
        return tool

    # Execute command and catch any exception
    def _execute(self, cmd):
        """
        Execute command
        """
        # Execute command, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            rc = -1

        # Return return code
        return rc
