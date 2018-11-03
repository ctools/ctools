#! /usr/bin/env python
# ==========================================================================
# [WHAT] script
#
# Copyright (C) [YEAR] [AUTHOR]
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
import sys
import gammalib
import ctools


# ================== #
# cscript_base class #
# ================== #
class cscript_base(ctools.cscript):
    """
    [WHAT] script
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor

        Parameters
        ----------
        argv : list of str
            List of IRAF command line parameter strings of the form
            ``parameter=3``.
        """
        # Initialise application by calling the base class constructor
        self._init_cscript(self.__class__.__name__, ctools.__version__, argv)

        # Return
        return

    # State methods for pickling
    def __getstate__(self):
        """
        Extend ctools.cscript __getstate__ method

        Returns
        -------
        state : dict
            Pickled instance
        """
        # Set pickled dictionary
        state = {'base' : ctools.cscript.__getstate__(self)}

        # Return pickled dictionary
        return state

    def __setstate__(self, state):
        """
        Extend ctools.cscript __setstate__ method

        Parameters
        ----------
        state : dict
            Pickled instance
        """
        # Set state
        ctools.cscript.__setstate__(self, state['base'])

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # TODO: Add code to query all relevant parameters

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return


    # Public methods
    def run(self):
        """
        Run the script
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()

        # TODO: Your code goes here

        # Return
        return

    def save(self):
        """
        Save something
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save something')

        # TODO: Your code goes here

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cscript_base(sys.argv)

    # Execute application
    app.execute()
