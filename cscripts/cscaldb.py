#! /usr/bin/env python
# ==========================================================================
# Show the content of the ctools calibration database
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
import sys
import glob
import os
import gammalib
import ctools


# ============= #
# cscaldb class #
# ============= #
class cscaldb(ctools.cscript):
    """
    Shows the content of the ctools calibration database.
    
    The ``cscaldb`` script writes the content of the calibration database
    into the ``cscaldb.log`` log file. If the ``debug`` parameter is set
    to ``yes`` the calibration database is also logged in the console.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        
        Parameters
        ----------
        argv : list of str
            List of IRAF command line parameter strings of the form
            ``parameter=3``.

        Raises
        ------
        TypeError
            An invalid number of command line arguments was provided.
        """
        # Set name
        self._name    = 'cscaldb'
        self._version = '1.1.0'

        # Initialise application by calling the appropriate class
        # constructor.
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    #
    # Query all user parameters
    def _get_parameters(self):

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    # Extract mission names from a calibration database
    def _get_missions(self, caldb):

        # Initialise mission list
        missions = []

        # Extract missions
        paths = glob.glob(caldb.rootdir()+'/data/*')
        paths.sort()
        for path in paths:
            missions.append(os.path.basename(path))

        # Sort missions
        missions.sort()

        # Return missions
        return missions

    # Extract instrument names from a calibration database
    def _get_instruments(self, caldb, mission):

        # Initialise instrument list
        instruments = []

        # Extract instruments
        paths       = glob.glob(caldb.rootdir()+'/data/'+mission+'/*')
        paths.sort()
        for path in paths:
            instruments.append(os.path.basename(path))

        # Sort instruments
        instruments.sort()

        # Return instruments
        return instruments

    # Extract response names from a calibrations FITS table and return them
    # in form of a list
    def _get_response_names(self, calibrations):

        # Initialise response name list
        names = []

        # Extract response names from calibrations and append them to the
        # response name list
        nrows = calibrations.length()
        ncols = calibrations.number()
        for row in range(nrows):
            for col in range(ncols):
                cal    = calibrations.string(row, col)
                istart = cal.find('NAME(')
                if istart != -1:
                    istop = cal.find(')')
                    name  = cal[5:istop]
                    if names.count(name) == 0:
                        names.append(name)
        
        # Sort response names
        names.sort()

        # Return response names
        return names


    # Public methods
    def run(self):
        """
        Run the script.
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()

        # Get the calibration database
        caldb = gammalib.GCaldb()

        # Extract mission names from the calibration database
        missions = self._get_missions(caldb)

        # Loop over missions
        for mission in missions:

            # Skip all non-CTA instruments
            if mission != 'cta':
                continue

            # Write mission into logger
            self._log_header1(gammalib.TERSE, 'Mission: '+mission)

            # Extract instruments
            instruments = self._get_instruments(caldb, mission)

            # Loop over instruments
            for instrument in instruments:

                # Write instrument into logger
                self._log_header3(gammalib.TERSE, 'Response functions in database "'+
                                                  instrument+'"')

                # Open calibration index file and retrieve calibrations
                filename = '/data/'+mission+'/'+instrument+'/caldb.indx'
                cifname  = caldb.rootdir() + filename
                fits     = gammalib.GFits(cifname)
                cif      = fits['CIF']
                caltable = cif['CAL_CBD']

                # Extract response names
                names = self._get_response_names(caltable)

                # Print response name
                if self._logTerse():
                    for name in names:
                        self._log(name+'\n')
                    self._log('\n')

        # Return
        return

    def execute(self):
        """
        Execute the script.
        """
        # Open logfile
        self.logFileOpen()

        # Run the script
        self.run()

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cscaldb(sys.argv)

    # Execute application
    app.execute()
