#! /usr/bin/env python
# ==========================================================================
# This script dumps all available calibrations into the console.
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
import gammalib
import ctools
import sys
import glob
import os


# ============= #
# cscaldb class #
# ============= #
class cscaldb(ctools.cscript):
    """
    The cscaldb class logs the content of the ctools calibration database
    into the log file. Optionally, if the "debug" parameter is set to "yes"
    the calibration database is also logged into the console.
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name    = "cscaldb"
        self._version = "1.0.0"

        # Initialise application
        if len(argv) == 0:
            ctools.cscript.__init__(self, self._name, self._version)
        elif len(argv) ==1:
            ctools.cscript.__init__(self, self._name, self._version, *argv)
        else:
            raise TypeError("Invalid number of arguments given.")

        # Set logger properties
        self.log_header()
        self.log.date(True)

        # Return
        return

    def __del__(self):
        """
        Destructor.
        """
        # Return
        return

    def _get_parameters(self):
        """
        Get parameters from parfile.
        """
        # Return
        return

    def _get_missions(self, caldb):
        """
        Extract mission names from a calibration database.

        Args:
            caldb: Calibration database.

        Returns:
            A list of mission names.
        """
        # Initialise mission list
        missions = []

        # Extract missions
        paths = glob.glob(caldb.rootdir()+"/data/*")
        paths.sort()
        for path in paths:
            missions.append(os.path.basename(path))

        # Sort missions
        missions.sort()

        # Return missions
        return missions

    def _get_instruments(self, caldb, mission):
        """
        Extract instrument names from a calibration database.

        Args:
            caldb:   Calibration database.
            mission: Mission name.

        Returns:
            A list of instrument names.
        """
        # Initialise instrument list
        instruments = []

        # Extract instruments
        paths       = glob.glob(caldb.rootdir()+"/data/"+mission+"/*")
        paths.sort()
        for path in paths:
            instruments.append(os.path.basename(path))

        # Sort instruments
        instruments.sort()

        # Return instruments
        return instruments

    def _get_response_names(self, calibrations):
        """
        Extract response names from a calibrations FITS table and return
        them in form of a list.

        Args:
            calibrations: FITS table containing the calibrations.

        Returns:
            A list of response names.
        """
        # Initialise response name list
        names = []

        # Extract response names from calibrations and append them to the
        # response name list
        nrows = calibrations.length()
        ncols = calibrations.number()
        for row in range(nrows):
            for col in range(ncols):
                cal    = calibrations.string(row, col)
                istart = cal.find("NAME(")
                if istart != -1:
                    istop = cal.find(")")
                    name  = cal[5:istop]
                    if names.count(name) == 0:
                        names.append(name)
        
        # Sort response names
        names.sort()

        # Return response names
        return names

    def execute(self):
        """
        Execute the script.
        """
        # Run the script
        self.run()

        # Return
        return

    def run(self):
        """
        Run the script.
        """
        # Switch screen logging on in debug mode
        if self.logDebug():
            self.log.cout(True)

        # Get parameters
        self._get_parameters()

        #  Write input parameters into logger
        if self.logTerse():
            self.log_parameters()
            self.log("\n")

        # Get the calibration database
        caldb = gammalib.GCaldb()

        # Extract mission names from the calibration database
        missions = self._get_missions(caldb)

        # Loop over missions
        for mission in missions:

            # Skip all non-CTA instruments
            if mission != "cta":
                continue

            # Write mission into logger
            if self.logTerse():
                self.log("\n")
                self.log.header1("Mission: "+mission)

            # Extract instruments
            instruments = self._get_instruments(caldb, mission)

            # Loop over instruments
            for instrument in instruments:

                # Write mission into logger
                if self.logTerse():
                    self.log.header3("Response functions in database \""+
                                     instrument+"\"")

                # Open calibration index file and retrieve calibrations
                filename = "/data/"+mission+"/"+instrument+"/caldb.indx"
                cifname  = caldb.rootdir() + filename
                fits     = gammalib.GFits(cifname)
                cif      = fits["CIF"]
                cals     = cif["CAL_CBD"]

                # Extract response names
                names = self._get_response_names(cals)

                # Print response name
                if self.logTerse():
                    for name in names:
                        self.log(name+"\n")
                    self.log("\n")

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Logs calibration database into the log file.
    """
    # Create instance of application
    app = cscaldb(sys.argv)

    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
