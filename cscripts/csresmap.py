#! /usr/bin/env python
# ==========================================================================
# Generates a residual map.
#
# Copyright (C) 2014-2016 Michael Mayer
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


# ============== #
# csresmap class #
# ============== #
class csresmap(ctools.cscript):
    """
    Generates a residual map.
    
    This class implements the creation of a residual map. It derives from
    the ctools.cscript class which provides support for parameter files,
    command line arguments, and logging. In that way the Python script
    behaves just as a regular ctool. 
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name
        self._name    = "csresmap"
        self._version = "1.1.0"

        # Initialise class members
        self._algorithm    = "SUB"
        self._resmap       = None
        self._modcube      = "NONE"
        self._use_maps     = False
        self._skip_binning = False
        self._publish      = False
        self._outfile      = ""
        self._use_maps     = False
        self._skip_binning = False
        self._edisp        = False
        self._xref         = 83.63
        self._yref         = 22.01
        self._emin         = 0.1
        self._emax         = 100.0
        self._enumbins     = 0
        self._ebinalg      = "LOG"
        self._coordsys     = "CEL"
        self._proj         = "CAR"
        self._nxpix        = 200
        self._nypix        = 200
        self._binsz        = 0.02
        self._publish      = False
        self._chatter      = 2
        self._clobber      = True
        self._debug        = False
        self._log_clients  = False

        # Initialise observation container from constructor arguments.
        self._obs, argv = self._set_input_obs(argv)

        # Initialise application by calling the appropriate class
        # constructor.
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Initialise some flags
        self._use_maps     = False
        self._skip_binning = False

        # First check if the inobs parameter is a counts cube
        if self._obs.size() == 0 and self["inobs"].filename() != "NONE":
            filename = gammalib.GFilename(self["inobs"].filename())
            if filename.is_fits():
                cta = gammalib.GCTAObservation()
                cta.load(filename)
                if cta.eventtype() == "CountsCube":
                    self._skip_binning = True

        # If we have a counts cube, then ask whether we also have a model
        if self._skip_binning:
            self._modcube = self["modcube"].filename()
            if self._modcube != "NONE":
                self._use_maps = True

        # If not two maps are given, proceed to set up observation
        if not self._use_maps:

            # Set observation if not done before
            if self._obs.size() == 0:
                self._require_inobs("csresmap.get_parameters()")
                self._obs = self._get_observations()

            # Check if we have exactly one binned CTA observation
            if self._obs.size() == 1:

                if self._obs[0].classname() == "GCTAObservation":
                    if self._obs[0].eventtype() == "CountsCube":                
                        # Skip ctbin step later on
                        self._skip_binning = True

            # Set models if we have none
            if self._obs.models().size() == 0:
                self._obs.models(self["inmodel"].filename())

            # Skip query for spatial parameters if a binning is provided in the observation
            if not self._skip_binning:
                # Read other parameters        
                self._xref     = self["xref"].real()
                self._yref     = self["yref"].real()
                self._emin     = self["emin"].real()
                self._emax     = self["emax"].real()
                self._enumbins = self["enumbins"].integer()
                self._ebinalg  = self["ebinalg"].string()
                self._coordsys = self["coordsys"].string()
                self._proj     = self["proj"].string()
                self._nxpix    = self["nxpix"].integer()
                self._nypix    = self["nypix"].integer()
                self._binsz    = self["binsz"].real()
                
        # Read energy dispersion flag
        self._edisp = self["edisp"].boolean()

        # Read algorithm
        self._algorithm = self["algorithm"].string()

        # Read standard parameters
        self._publish = self["publish"].boolean()
        self._chatter = self["chatter"].integer()
        self._clobber = self["clobber"].boolean()
        self._debug   = self["debug"].boolean()

        # Read ahead output parameters
        if (self._read_ahead()):
            self["outmap"].filename()

        # Write input parameters into logger
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

        # Write observation into logger
        if self._logTerse():
            self._log("\n")
            self._log.header1("Observation")
            self._log(str(self._obs))
            self._log("\n")

        # Use input file directly if given
        if self._use_maps:
            countmap = gammalib.GSkyMap(self["inobs"].filename())
            modelmap = gammalib.GSkyMap(self._modcube)

        else:

            # ...
            if self._skip_binning:
                cta_counts_cube = gammalib.GCTAEventCube(self._obs[0].events().clone())   

            # ...
            else:

                # Write header
                if self._logTerse():
                    self._log("\n")
                    self._log.header1("Generate binned map (ctbin)")

                # Create countsmap
                bin = ctools.ctbin(self._obs)
                bin["nxpix"].integer(self._nxpix)
                bin["nypix"].integer(self._nypix)
                bin["proj"].string(self._proj)
                bin["coordsys"].string(self._coordsys)
                bin["xref"].real(self._xref)
                bin["yref"].real(self._yref)
                bin["enumbins"].integer(self._enumbins)
                bin["ebinalg"].string(self._ebinalg)
                bin["emin"].real(self._emin)
                bin["emax"].real(self._emax)
                bin["binsz"].real(self._binsz)
                bin["chatter"].integer(self._chatter)
                bin["clobber"].boolean(self._clobber)
                bin["debug"].boolean(self._debug)
                bin.run()

                # Retrieve counts cube
                cta_counts_cube = bin.cube()

            # Assign GCTAEventCube to skymap
            countmap = cta_counts_cube.counts()

            # Write header
            if self._logTerse():
                self._log("\n")
                self._log.header1("Generate model map (ctmodel)")

            # Create model map
            model = ctools.ctmodel(self._obs)
            model.cube(cta_counts_cube)
            model["chatter"].integer(self._chatter)
            model["clobber"].boolean(self._clobber)
            model["debug"].boolean(self._debug)
            model["edisp"].boolean(self._edisp)
            model.run()

            # Get model map into GSkyMap object
            modelmap = model.cube().counts().copy()

        # Store counts map as residual map. Note that we need a
        # special construct here to avoid memory leaks. This seems
        # to be a SWIG feature as SWIG creates a new object when
        # calling bin.cube()
        #residualmap = bin.cube().counts()
        self._resmap = countmap.copy()
        self._resmap.stack_maps()
        modelmap.stack_maps()

        # Continue calculations depending on given algorithm
        if self._algorithm == "SUB":

            # Subtract maps 
            self._resmap -= modelmap

        elif self._algorithm == "SUBDIV":

            # Subtract and divide by model map
            self._resmap -= modelmap
            self._resmap /= modelmap
            #for pixel in modelmap:
            #    if pixel != 0.0:
            #        pixel = 1.0/pixel
            #self._resmap *= modelmap

        elif self._algorithm == "SUBDIVSQRT":

            # subtract and divide by sqrt of model map
            self._resmap -= modelmap
            self._resmap /= modelmap.sqrt()
            #for pixel in modelmap:
            #    if pixel != 0.0:
            #        pixel = 1.0/math.sqrt(pixel)
            #self._resmap *= modelmap

        else:

            # Raise error if algorithm is unkown
            raise TypeError("Algorithm \""+self._algorithm+"\" not known")

        # Optionally publish map
        if self._publish:
            self.publish()

        # Return
        return

    def execute(self):
        """
        Execute the script.
        """
        # Open logfile
        self.logFileOpen()

        # Read ahead output parameters
        self._read_ahead(True)

        # Run the script
        self.run()

        # Save residual map
        self.save()

        # Return
        return

    def save(self):
        """
        Save residual map.
        """
        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Save residual map")

        # Get outmap parameter
        outmap = self["outmap"].filename()
        
        # Continue only filename and residual map are valid
        if self._resmap != None:

            # Log file name
            if self._logTerse():
                self._log(gammalib.parformat("Residual map file"))
                self._log(outmap.url())
                self._log("\n")

            # Save residual map
            self._resmap.save(outmap, self._clobber)

        # Return
        return

    def publish(self, name=""):
        """
        Publish residual map.

        Kwargs:
            name: Name of residual map.
        """
        # Write header
        if self._logTerse():
            self._log("\n")
            self._log.header1("Publish residual map")

        # Continue only if residual map is valid
        if self._resmap != None:
        
            # Set default name is user name is empty
            if not name:
                user_name = self.name
            else:
                user_name = name

            # Log map name
            if self._logTerse():
                self._log("Publish residual map \""+user_name+"\".\n")

            # Publish map
            self._resmap.publish(user_name)

        # Return
        return

    def models(self, models):
        """
        Set model.
        """
        # Copy models
        self._obs.models(models.clone())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csresmap(sys.argv)

    # Execute application
    app.execute()
