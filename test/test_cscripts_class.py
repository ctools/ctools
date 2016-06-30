#! /usr/bin/env python
# ==========================================================================
# cscripts test base class
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


# ======================= #
# Test class for cscripts #
# ======================= #
class cscripts_test(gammalib.GPythonTestSuite):
    """
    Test class for scripts
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Set script name
    def _script(self, name):
        """
        Set script name
        """
        # Kluge to set the command. The installed version has no README file
        if os.path.isfile('README.md'):
            script = '../cscripts/'+name+'.py'
        else:
            script = name

        # Return script name
        return script

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
            pass

        # Return return code
        return rc
