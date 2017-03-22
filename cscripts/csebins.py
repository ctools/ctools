#! /usr/bin/env python
# ==========================================================================
# Generates energy boundaries for stacked analysis
#
# Copyright (C) 2017 Juergen Knoedlseder
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
import math
import gammalib
import ctools
#from cscripts import obsutils


# ============= #
# csebins class #
# ============= #
class csebins(ctools.cscript):
    """
    Generates energy boundaries for stacked analysis
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Set name and version
        self._name    = 'csebins'
        self._version = '1.0.0'

        # Set members
        self._ebounds = gammalib.GEbounds()

        # Initialise observation container from constructor arguments.
        self._obs, argv = self._set_input_obs(argv)

        # Initialise script by calling the appropriate class constructor.
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Initialise observation container if it is currently empty
        if self._obs.size() == 0:

            # Query "inobs" parameter
            filename = self['inobs'].filename()

            # If an observation definition file was provided then load it,
            # otherwise build a single observation with response information
            if self._is_valid_filename(filename):
                self._obs = gammalib.GObservations(filename)
            else:
                cta   = gammalib.GCTAObservation()
                caldb = gammalib.GCaldb('cta', self['caldb'].string())
                rsp   = gammalib.GCTAResponseIrf(self['irf'].string(), caldb)
                cta.response(rsp)
                self._obs.append(cta)

        # Query input parameters
        self['emin'].real()
        self['emax'].real()
        self['aeffthres'].real()
        self['bkgthres'].real()

        # Query ahead output model filename
        if self._read_ahead():
            self['outfile'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _insert_energy(self, responses, energy, comment):
        """
        Set energy boundaries for one CTA observation

        Parameters
        ----------
        responses : list of dict
            List of response dictionaries
        energy : float
            Energy to insert (TeV)
        comment : str
            Reason for energy insertion
        """
        # Set thresholds of all response components
        logE = math.log10(energy)
        for rsp in responses:
            rsp['thres'] = rsp['irf'](logE, 0.0, 0.0)

        # Convert energy into gammalib energy
        eng = gammalib.GEnergy(energy, 'TeV')
  
        # Set insertion message
        msg = str(energy)+' TeV'
        if len(comment) > 0:
            msg += ' ('+comment+')'
  
        # If energy boundaries are empty then append an energy interval
        # with zero size ...
        if self._ebounds.size() == 0:

            # Append energy interval with zero size
            self._ebounds.append(eng, eng)

            # Log first energy boundary
            self._log_value(gammalib.NORMAL, 'First boundary', msg)

        # ... otherwise append an interval
        elif self._ebounds.size() > 0:

            # If current boundaries have zero size then recover boundary energy,
            # remove boundary, and append a boundary of non-zero size
            if self._ebounds.emin() == self._ebounds.emax():

                # Recover old energy, remove boundary and append new boundary
                emin = self._ebounds.emin().copy() # Because emin becomes zero after remove
                emax = eng
                self._ebounds.remove(0)
                self._ebounds.append(emin, emax)

                # Log other energy boundary
                self._log_value(gammalib.NORMAL, 'Boundary', msg)

            else:

                # Use old emax as new emin. Only append if interval is non zero
                emin = self._ebounds.emax()
                if emin < eng:
                    self._ebounds.append(emin, eng)

                    # Log other energy boundary
                    self._log_value(gammalib.NORMAL, 'Boundary', msg)

        # Return
        return

    def _set_ebounds(self, obs):
        """
        Set energy boundaries for one CTA observation

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            CTA observation
        """
        # Get parameters
        emin      = self['emin'].real()
        emax      = self['emax'].real()
        aeffthres = self['aeffthres'].real()
        bkgthres  = self['bkgthres'].real()

        # Build list of response dictionaries
        responses = []
        
        # Loop over all observations
        for run in obs:

            # Get response name
            mission    = run.response().caldb().mission()
            instrument = run.response().caldb().instrument()
            rspname    = run.response().rspname()
            name       = '%s::%s::%s' % (mission, instrument, rspname)

            # If response name exists already then examine next observation
            exists = False
            for rsp in responses:
                if rsp['name'] == name:
                    exists = True
                    break
            if exists:
                continue

            # Append effective area and background template to list
            aeff = {'name': name, 'type': 'Aeff',
                    'irf': run.response().aeff(), 'thres': 0.0}
            bkg  = {'name': name, 'type': 'Background',
                    'irf': run.response().background(), 'thres': 0.0}
            responses.append(aeff)
            responses.append(bkg)

            # Log response name
            self._log_value(gammalib.NORMAL, 'Append response', '%s (%s) [%s]' %
                            (rspname, instrument, mission))

        # Setup energy vector in log10(TeV)
        logEs = [math.log10(1.0e-3*float(i)) for i in range(int(emin*1000),int(emax*1000))]

        # Loop over all energies
        for logE in logEs:

            # Loop over all response components
            for rsp in responses:

                # Get IRF value
                irf = rsp['irf'](logE, 0.0, 0.0)

                # Get threshold
                if rsp['type'] == 'Aeff':
                    threshold = aeffthres
                else:
                    threshold = bkgthres

                # If threshold is zero and response is positive then insert
                # an energy ...
                if rsp['thres'] == 0.0 and irf > 0.0:
                    self._insert_energy(responses, math.pow(10.0, logE), rsp['type'])
                    break

                # ... otherwise if threshold is not zero then check whether
                # the fractional change exceeds the threshold
                elif rsp['thres'] != 0.0:
                    f = irf/rsp['thres']
                    if 1.0-f > threshold or f-1.0 > threshold:
                        self._insert_energy(responses, math.pow(10.0, logE), rsp['type'])
                        break

        # Insert energy at end of range
        self._insert_energy(responses, emax, 'End of range')

        # Log number of energy boundaries
        self._log_value(gammalib.NORMAL, 'Number of boundaries', self._ebounds.size())

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

        # Write observation into logger
        self._log_observations(gammalib.NORMAL, self._obs, 'Input observation')

        # Write header
        self._log_header1(gammalib.TERSE, 'Define energy boundaries')

        # Define energy boundaries for observation container
        self._set_ebounds(self._obs)

        # Return
        return

    def execute(self):
        """
        Execute the script
        """
        # Open logfile
        self.logFileOpen()

        # Read ahead output parameters
        self._read_ahead(True)

        # Run the scrip
        self.run()

        # Save energy boundaries
        self.save()

        # Return
        return

    def save(self):
        """
        Save energy boundaries
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save energy boundaries')

        # Get energy boundaries filename
        outfile = self['outfile'].filename()

        # Log file name
        self._log_value(gammalib.NORMAL, 'Energy boundaries file', outfile.url())

        # Save energy boundaries
        self._ebounds.save(outfile, self._clobber())

        # Return
        return

    def ebounds(self):
        """
        Return energy boundaries

        Returns
        -------
        ebounds : `~gammalib.GEbounds()`
            Energy boundaries
        """
        # Return
        return self._ebounds


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csebins(sys.argv)

    # Execute application
    app.execute()
