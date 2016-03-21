#! /usr/bin/env python
# ==========================================================================
# Generates model cube for Fermi/LAT data.
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
import sys


# ============= #
# lsmodel class #
# ============= #
class lsmodel(gammalib.GApplication):
    """
    Generates model cube for Fermi/LAT data.
    
    This class implements the model map generation script. It derives from
    the GammaLib::GApplication class which provides support for parameter
    files, command line arguments, and logging. In that way the Python
    script behaves just as a regular ctool.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name and version
        self._name    = "lsmodel"
        self._version = "1.1.0"

        # Initialise application
        if len(argv) == 0:
            gammalib.GApplication.__init__(self, self._name, self._version)
        elif len(argv) ==1:
            gammalib.GApplication.__init__(self, self._name, self._version, *argv)
        else:
            raise TypeError("Invalid number of arguments given.")

        # Set logger properties
        self._log_header()
        self._log.date(True)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile.
        """
        # Query user parameters
        self["cntmap"].filename()
        self["expmap"].filename()
        self["ltcube"].filename()
        self["srcmdl"].filename()
        self["caldb"].string()
        self["irf"].string()

        # Read ahead output parameters
        if (self._read_ahead()):
            self["outfile"].filename()

        #  Write input parameters into logger
        if self._logTerse():
            self._log_parameters()
            self._log("\n")

        # Return
        return


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

        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Generate model cube")

        # Get user parameters
        cntmap  = self["cntmap"].filename()
        expmap  = self["expmap"].filename()
        ltcube  = self["ltcube"].filename()
        srcmdl  = self["srcmdl"].filename()
        caldb   = self["caldb"].string()
        irf     = self["irf"].string()

        # Set LAT observation
        obs = gammalib.GLATObservation()
        obs.load_binned(cntmap, expmap, ltcube)
        obs.response(irf, caldb)

        # Set model (use first model for the moment)
        models = gammalib.GModels(srcmdl)
        model  = models[0]

        # Get deep copy of GLATEventCube from observation
        cube = gammalib.GLATEventCube(obs.events())

        # Fill cube with model values
        for i in range(cube.size()):
            event = cube[i]
            value = model.eval(event, obs) * event.size()
            event.counts(value)
            cube[i] = event

        # Log cube
        if self._logTerse():
            self._log(str(cube))
            self._log("\n")

        # Return
        return

    def execute(self):
        """
        Execute the script.
        """
        # Open logfile
        self._logFileOpen()

        # Read ahead output parameters
        self._read_ahead(True)

        # Run the script
        self.run()

        # Save event cube
        self.save()

        # Return
        return

    def save(self):
        """
        Save model cube.
        """
        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Save model cube")

        # Get outfile parameter
        outfile = self["outfile"].filename()
        
        # Log file name
        if self._logTerse():
            self._log(gammalib.parformat("Model cube"))
            self._log(outfile.url())
            self._log("\n")

        # Save event cube
        self._cube.save(outfile, self._clobber())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = lsmodel(sys.argv)

    # Execute application
    app.execute()
