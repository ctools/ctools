#! /usr/bin/env python
# ==========================================================================
# Generate pulsar phase bins
#
# Copyright (C) 2022 Juergen Knoedlseder
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


# =============== #
# compulbin class #
# =============== #
class compulbin(ctools.csobservation):
    """
    Generate pulsar phase bins
    """
    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the base class constructor
        self._init_csobservation(self.__class__.__name__, ctools.__version__, argv)

        # Initialise members
        self._select   = gammalib.GCOMSelection()
        self._pnumbins = 0
        self._phasebin = 1.0
        self._bins     = []
        self._fits     = None

        # Return
        return

    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Set observation if not done before
        if self.obs().is_empty():
            self.obs().load(self['inobs'].filename())

        # Query parameters
        self['psrname'].string()
        self['ephemerides'].filename()
        self['emin'].real()
        self['emax'].real()
        self._pnumbins = self['pnumbins'].integer()
        self['armmin'].real()
        self['armmax'].real()
        self['phimin'].real()
        self['phimax'].real()
        self['zetamin'].real()
        self['fpmtflag'].integer()
        self['psdmin'].integer()
        self['psdmax'].integer()

        # Get D1 and D2 module usage strings
        d1use = self['d1use'].string()
        d2use = self['d2use'].string()

        # Check D1 and D2 module usage strings
        if len(d1use) != 7:
            msg = 'Incorrect length %d of D1 usage string. String needs to have 7 digits.' % \
                  len(d1use)
            raise RuntimeError(msg)
        if len(d2use) != 14:
            msg = 'Incorrect length %d of D2 usage string.  String needs to have 14 digits.' % \
                  len(d2use)
            raise RuntimeError(msg)

        # Query ahead output model filename
        if self._read_ahead():
            self['outfile'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _set_selection_set(self):
        """
        Set selection set
        """
        # Write header
        self._log_header1(gammalib.NORMAL, 'Set selection set')

        # Clear selection set
        self._select.clear()

        # Get pulsar ephemerides
        pulsar = gammalib.GPulsar(self['ephemerides'].filename(), self['psrname'].string())

        # Set pulsar ephemerides
        self._select.pulsar(pulsar)

        # Set PSD interval
        self._select.psd_min(self['psdmin'].integer())
        self._select.psd_max(self['psdmax'].integer())

        # Set handling of D2 modules with failed PMT flag
        self._select.fpmtflag(self['fpmtflag'].integer())

        # Set D1 module usage flags
        d1use = self['d1use'].string()
        for i in range(7):
            if d1use[i] == '1':
                self._select.use_d1(i,True)
            else:
                self._select.use_d1(i,False)

        # Set D2 module usage flags
        d2use = self['d2use'].string()
        for i in range(14):
            if d2use[i] == '1':
                self._select.use_d2(i,True)
            else:
                self._select.use_d2(i,False)

        # Log selection set
        self._log_string(gammalib.NORMAL, str(self._select))

        # Return
        return

    def _bin_observation(self, obs, zetamin=5.0):
        """
        Bin events in phase bins

        Parameters
        ----------
        obs : `~gammalib.GCOMObservation`
            Unbinned COMPTEL observation
        zetamin : float, optional
            Minimum zeta angle
        """
        # Initialise event index
        ievent = 0

        # Initialise statistics
        num_used_superpackets    = 0
        num_skipped_superpackets = 0
        num_processed            = 0
        num_event_outside_sp     = 0
        num_energy_too_low       = 0
        num_energy_too_high      = 0
        num_phibar_too_small     = 0
        num_phibar_too_large     = 0
        num_eha_too_small        = 0
        num_arm_too_small        = 0
        num_arm_too_large        = 0
        num_bin_too_small        = 0
        num_bin_too_large        = 0
        num_used_events          = 0

        # Get selection parameters
        emin   = gammalib.GEnergy(self['emin'].real(), 'MeV')
        emax   = gammalib.GEnergy(self['emax'].real(), 'MeV')
        armmin = self['armmin'].real()
        armmax = self['armmax'].real()
        phimin = self['phimin'].real()
        phimax = self['phimax'].real()

        # Get Good Time Intervals and reduce them to the validity range
        # of the ephemerides
        tim = obs.tim()
        tim.reduce(self._select.pulsar().validity())

        # Loop over Orbit Aspect Data
        for oad in obs.oads():

            # Skip superpacket if it is not fully enclosed within the COMPTEL
            # Good Time Interval
            if not tim.contains(oad.tstart()) or not tim.contains(oad.tstop()):
                num_skipped_superpackets += 1
                continue

            # Increment superpacket counter
            num_used_superpackets += 1

            # Collect all events in the superpacket
            while ievent < obs.events().size():

                # Get event
                event = obs.events()[ievent]

                # Break event loop if end of superpacket was reached
                if event.time() > oad.tstop():
                    break

                # Increment event counter
                ievent += 1

                # Skip event if it lies before the superpacket start
                if event.time() < oad.tstart():
                    num_event_outside_sp += 1
                    continue

                # Increment number of processed events
                num_processed += 1

                # Skip event if it lies outside energy range
                if event.energy() < emin:
                    num_energy_too_low += 1
                    continue
                elif event.energy() > emax:
                    num_energy_too_high += 1
                    continue

                # Apply event selection
                if not self._select.use_event(event):
                    continue

                # Skip event if it lies outside Phibar range
                if event.phibar() < phimin:
                    num_phibar_too_small += 1
                    continue
                elif event.phibar() > phimax:
                    num_phibar_too_large += 1
                    continue

                # Skip event if it comes from the Earth horizon
                ehamin = event.phibar() + zetamin
                if event.eha() < ehamin:
                    num_eha_too_small += 1
                    continue

                # Get pulsar ephemeris for event time
                ephemeris = self._select.pulsar().ephemeris(event.time())

                # Compute Phigeo
                phigeo = event.dir().dir().dist_deg(ephemeris.dir())

                # Make ARM selection
                arm = event.phibar() - phigeo
                if arm < armmin:
                    num_arm_too_small += 1
                    continue
                elif arm > armmax:
                    num_arm_too_large += 1
                    continue

                # Convert event time to Solar System Barycentre time
                tdelta = obs.bvcs().tdelta(ephemeris.dir(), event.time())
                time   = event.time() + tdelta

                # Compute pulsar phase
                phase = ephemeris.phase(time)

                # Fill event in phase bin
                ibin = int(phase / self._phasebin)
                if ibin < 0:
                    num_bin_too_small += 1
                elif ibin >= self._pnumbins:
                    num_bin_too_large += 1
                else:
                    self._bins[ibin] += 1.0
                    num_used_events  += 1

            # Break superpacket loop if there are no more events
            if ievent >= obs.events().size():
                break

        # If not all events are exhausted then add them to the events outside
        # superpackets
        if ievent < obs.events().size():
            num_event_outside_sp += (obs.events().size() - ievent)

        # Log statistics
        self._log_header3(gammalib.NORMAL, 'Event binning statistics')
        self._log_value(gammalib.NORMAL, 'Number of superpackets', obs.oads().size())
        self._log_value(gammalib.NORMAL, 'Used superpackets', num_used_superpackets)
        self._log_value(gammalib.NORMAL, 'Skipped superpackets', num_skipped_superpackets)
        self._log_value(gammalib.NORMAL, 'Total number of events', obs.events().size())
        self._log_value(gammalib.NORMAL, 'Processed events', num_processed)
        self._log_value(gammalib.NORMAL, 'Used events', num_used_events)
        self._log_value(gammalib.NORMAL, 'Events outside superpacket', num_event_outside_sp)
        self._log_value(gammalib.NORMAL, 'Energy too low', num_energy_too_low)
        self._log_value(gammalib.NORMAL, 'Energy too high', num_energy_too_high)
        self._log_value(gammalib.NORMAL, 'Phibar too small', num_phibar_too_small)
        self._log_value(gammalib.NORMAL, 'Phibar too large', num_phibar_too_large)
        self._log_value(gammalib.NORMAL, 'Earth hor. angle too small', num_eha_too_small)
        self._log_value(gammalib.NORMAL, 'ARM too small', num_arm_too_small)
        self._log_value(gammalib.NORMAL, 'ARM too large', num_arm_too_large)
        self._log_value(gammalib.NORMAL, 'Phase bin index too small', num_bin_too_small)
        self._log_value(gammalib.NORMAL, 'Phase bin index too large', num_bin_too_large)

        # Return
        return

    def _create_fits(self):
        """
        Creates FITS file with phase bins
        """
        # Create FITS table columns
        phase_min = gammalib.GFitsTableFloatCol('PHASE_MIN', self._pnumbins)
        phase_max = gammalib.GFitsTableFloatCol('PHASE_MAX', self._pnumbins)
        nevents   = gammalib.GFitsTableDoubleCol('NEVENTS', self._pnumbins)
        nevents.unit('counts')

        # Fill FITS table columns
        for i in range(self._pnumbins):
            phase_min[i] = i * self._phasebin
            phase_max[i] = (i+1) * self._phasebin
            nevents[i]   = self._bins[i]

        # Create FITS Table with extension "PHASECURVE"
        table = gammalib.GFitsBinTable(self._pnumbins)
        table.extname('PHASECURVE')

        # Add keywords
        table.card('INSTRUME', 'CGRO', 'Name of Instrument')
        table.card('TELESCOP', 'COMPTEL', 'Name of Telescope')

        # Stamp header
        self._stamp(table)

        # Add script keywords
        table.card('PSRNAME',  self['psrname'].string(), 'Pulsar name')
        table.card('EPHEM',    self['ephemerides'].filename().url(), 'Ephemerides file')
        table.card('EMIN',     self['emin'].real(), '[MeV] Minimum energy')
        table.card('EMAX',     self['emax'].real(), '[MeV] Maximum energy')
        table.card('ARMMIN',   self['armmin'].real(), '[deg] Minimum angular resolution measure')
        table.card('ARMMAX',   self['armmax'].real(), '[deg] Maximum angular resolution measure')
        table.card('PHIMIN',   self['phimin'].real(), '[deg] Minimum Phibar')
        table.card('PHIMAX',   self['phimax'].real(), '[deg] Maximum Phibar')
        table.card('ZETAMIN',  self['zetamin'].real(), '[deg] Minimum zeta angle')
        table.card('FPMTFLAG', self['fpmtflag'].integer(), 'Handling of failed PMTs')
        table.card('PSDMIN',   self['psdmin'].integer(), 'Minimum PSD channel')
        table.card('PSDMAX',   self['psdmax'].integer(), 'Maximum PSD channel')
        table.card('D1USE',   self['d1use'].string(), 'D1 modules usage')
        table.card('D2USE',   self['d2use'].string(), 'D2 modules usage')

        # Append filled columns to fits table
        table.append(phase_min)
        table.append(phase_max)
        table.append(nevents)

        # Create the FITS file now
        self._fits = gammalib.GFits()
        self._fits.append(table)

        # Return table
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

        # Initialise phase bins
        self._phasebin = 1.0 / self._pnumbins
        self._bins     = [0.0 for i in range(self._pnumbins)]

        # Log header
        self._log_header1(gammalib.NORMAL, 'Input observations')

        # Log input observations
        self._log_string(gammalib.NORMAL, str(self.obs()))

        # Set selection set
        self._set_selection_set()

        # Write header
        self._log_header1(gammalib.NORMAL, 'Generate pulsar phase bins')

        # Loop over all input observations
        for obs in self.obs():

            # Write header
            self._log_header2(gammalib.NORMAL, self._get_obs_header(obs))

            # Log observation
            self._log_string(gammalib.NORMAL, str(obs))

            # Generate phase bins
            self._bin_observation(obs)

        # Create FITS file
        self._create_fits()

        # Return
        return

    def save(self):
        """ 
        Save phase bin file
        """
        # Write header
        self._log_header1(gammalib.TERSE, 'Save phase bin file')

        # Continue only if FITS file is valid
        if self._fits != None:

            # Get output filename
            outfile = self['outfile'].filename()

            # If file exists and clobber flag is false then raise an exception
            if outfile.exists() and not self['clobber'].boolean():
                msg = ('Cannot save "'+outfile.url()+'": File already exists. '
                       'Use parameter clobber=yes to allow overwriting of files.')
                raise RuntimeError(msg)

            # ... otherwise log filename and save file
            else:
                # Log filename
                self._log_value(gammalib.NORMAL, 'Phase bin file', outfile.url())

                # Save phase bins
                self._fits.saveto(outfile, self['clobber'].boolean())

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = compulbin(sys.argv)

    # Execute application
    app.execute()
