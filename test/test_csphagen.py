#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csphagen script.
#
# Copyright (C) 2017 Luigi Tibaldo
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
import cscripts
from testing import test


# ================================ #
# Test class for csphagen script #
# ================================ #
class Test(test):
    """
    Test class for csphagen script

    This test class makes unit tests for the csphagen script by using it
    from the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)
        self._myevents1 = self._datadir + '/crab_offaxis1.fits'
        self._myevents2 = self._datadir + '/crab_offaxis2.fits'
        self._exclusion = self._datadir + '/crab_exclusion.fits'
        self._nreg_with_excl = 5
        self._nreg_wo_excl = 8
        self._nreg_mul = [self._nreg_with_excl, 6]
        # number of expected background regions with/wo exclusion,
        # and for two different runs

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('csphagen')

        # Append tests
        self.append(self._test_cmd, 'Test csphagen on command line')
        self.append(self._test_python, 'Test csphagen from Python')

        # Return
        return

    # Test csphagen on command line
    def _test_cmd(self):
        """
        Test csphagen on the command line.
        """
        # Set script name
        csphagen = self._script('csphagen')
        nbins = 120

        cmd = csphagen + 'inobs="' + self._myevents1 + \
              '" caldb="' + self._caldb + '" irf="' + self._irf + \
              '" ebinalg=LOG emin=0.1 emax=100. enumbins="' + nbins + \
              '" coordsys=CEL' + ' ra=83.633 dec=22.0145 rad=0.2 stack=no exclusion="' + \
              self._exclusion + '" outroot=genpha_cmd1' + \
              ' logfile="csphagen_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
                         'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
                         'Check successful execution from command line')

        # Check output files
        self._check_output('genpha_cmd1_', nbins, self._nreg_with_excl)
        self._check_outobs('genpha_cmd1', 1)

        cmd = csphagen + 'inobs="events_that_do_not_exist.fits"' + \
              '" caldb="' + self._caldb + '" irf="' + self._irf + \
              '" ebinalg=LOG emin=0.1 emax=100. enumbins="' + nbins + \
              '" coordsys=CEL' + ' ra=83.633 dec=22.0145 rad=0.2 stack=no exclusion="' + \
              self._exclusion + '" outroot=genpha_cmd2' + \
              ' logfile="csphagen_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
                         'Check invalid input file when executed from command line')

        # Check cslightcrv --help
        self._check_help(csphagen)

        # Return
        return

    def _test_python(self):
        """
        Test csphagen from Python
        """
        nbins = 120

        # Same test as from command line
        phagen = cscripts.csphagen()
        phagen['inobs'] = self._myevents1
        phagen['caldb'] = self._caldb
        phagen['irf'] = self._irf
        phagen['ebinalg'] = 'LOG'
        phagen['emin'] = 0.1
        phagen['emax'] = 100.
        phagen['enumbins'] = nbins
        phagen['coordsys'] = 'CEL'
        phagen['ra'] = 83.633
        phagen['dec'] = 22.0145
        phagen['rad'] = 0.2
        phagen['stack'] = False
        phagen['exclusion'] = self._exclusion
        phagen['outroot'] = 'genpha_py1'
        phagen['logfile'] = 'csphagen_py1.log'
        phagen['chatter'] = 1

        # Run script
        phagen.execute()

        # Check outout
        self._check_output('genpha_py1_', nbins, self._nreg_with_excl)
        self._check_outobs('genpha_py1', 1)

        # Now test without exclusion region
        phagen = cscripts.csphagen()
        phagen['inobs'] = self._myevents1
        phagen['caldb'] = self._caldb
        phagen['irf'] = self._irf
        phagen['ebinalg'] = 'LOG'
        phagen['emin'] = 0.1
        phagen['emax'] = 100.
        phagen['enumbins'] = nbins
        phagen['coordsys'] = 'CEL'
        phagen['ra'] = 83.633
        phagen['dec'] = 22.0145
        phagen['rad'] = 0.2
        phagen['stack'] = False
        phagen['outroot'] = 'genpha_py2'
        phagen['logfile'] = 'csphagen_py2.log'
        phagen['chatter'] = 1

        # Run script
        phagen.execute()

        # Check outout
        self._check_output('genpha_py2_', nbins, self._nreg_wo_excl)
        self._check_outobs('genpha_py2', 1)

        # Test with multiple input observations, no stacking
        # Create observation container
        obs = gammalib.GObservations()
        for s, events in enumerate([self._myevents1, self._myevents2]):
            run = gammalib.GCTAObservation(events)
            run.id(str(s + 1))
            run.response(self.irf, self.caldb)
            obs.append(run)

        # Setup csphagen
        phagen = cscripts.csphagen(obs)
        phagen['ebinalg'] = 'LOG'
        phagen['emin'] = 0.1
        phagen['emax'] = 100.
        phagen['enumbins'] = nbins
        phagen['coordsys'] = 'CEL'
        phagen['ra'] = 83.633
        phagen['dec'] = 22.0145
        phagen['rad'] = 0.2
        phagen['stack'] = False
        phagen['exclusion'] = self._exclusion
        phagen['outroot'] = 'genpha_py3'
        phagen['logfile'] = 'csphagen_py3.log'
        phagen['chatter'] = 1

        # Run script
        phagen.execute()

        # Check outout
        for s in range(2):
            self._check_output('genpha_py3_' + str(s + 1), nbins,
                               self._nreg_mul[s])
        self._check_outobs('genpha_py3', 2)

        return

    def _check_ebounds(self, table, bins):
        """
        Check EBOUNDS table
        """
        cols = ['E_MIN', 'E_MAX']

        self.test_value(table.ncols(), len(cols),
                        'Check for %d columns in light curve FITS table' % len(
                            cols))
        self.test_value(table.nrows(), bins,
                        'Check for %d rows in light curve FITS table' % bins)
        for col in cols:
            self.test_assert(table.contains(col),
                             'FITS file contains "' + col + '" column')

        return

    def _check_pha(self, filename, bins):
        """
        Check PHA file
        """
        # Expected column names
        cols = ['CHANNEL', 'COUNTS', 'STAT_ERR', 'SYS_ERR', 'QUALITY',
                'GROUPING', 'AREASCAL', 'BACKSCAL']

        # Open FITS file
        fits = gammalib.GFits(filename)

        # Check FITS file structure
        self.test_value(fits.size(), 3,
                        'Check for 3 extensions in light curve FITS file')
        self.test_assert(fits.contains('SPECTRUM'),
                         'FITS file contains "SPECTRUM" extension')
        self.test_assert(fits.contains('EBOUNDS'),
                         'FITS file contains "EBOUNDS" extension')

        # Get SPECTRUM table
        table = fits['SPECTRUM']

        # Check FITS table structure
        self.test_value(table.ncols(), len(cols),
                        'Check for %d columns in light curve FITS table' % len(
                            cols))
        self.test_value(table.nrows(), bins,
                        'Check for %d rows in light curve FITS table' % bins)
        for col in cols:
            self.test_assert(table.contains(col),
                             'FITS file contains "' + col + '" column')

        # Check SPECTRUM table
        table = fits['EBOUNDS']
        self._check_ebounds(table, nbins)

        # Close FITS file
        fits.close()

        # Return
        return

    def _check_arf(self, filename, bins):
        """
        Check ARF file
        """
        # Expected column names
        cols = ['ENERG_LO', 'ENERG_HI', 'SPECRESP', 'BACKGROUND', 'ALPHA']

        # Open FITS file
        fits = gammalib.GFits(filename)

        # Check FITS file structure
        self.test_value(fits.size(), 2,
                        'Check for 2 extensions in light curve FITS file')
        self.test_assert(fits.contains('SPECRESP'),
                         'FITS file contains "SPECRESP" extension')

        # Get SPECRESP table
        table = fits['SPECRESP']

        # Check FITS table structure
        self.test_value(table.ncols(), len(cols),
                        'Check for %d columns in light curve FITS table' % len(
                            cols))
        self.test_value(table.nrows(), bins,
                        'Check for %d rows in light curve FITS table' % bins)
        for col in cols:
            self.test_assert(table.contains(col),
                             'FITS file contains "' + col + '" column')

        # Close FITS file
        fits.close()

        # Return
        return

    def _check_rmf(self, filename, bins):
        """
        Check RMF file
        """
        # Expected column names
        cols = ['ENERG_LO', 'ENERG_HI', 'N_GRP', 'F_CHAN', 'N_CHAN', 'MATRIX']

        # Open FITS file
        fits = gammalib.GFits(filename)

        # Check FITS file structure
        self.test_value(fits.size(), 3,
                        'Check for 3 extensions in light curve FITS file')
        self.test_assert(fits.contains('MATRIX'),
                         'FITS file contains "MATRIX" extension')
        self.test_assert(fits.contains('EBOUNDS'),
                         'FITS file contains "EBOUNDS" extension')

        # Get MATRIX table
        table = fits['MATRIX']

        # Check FITS table structure
        self.test_value(table.ncols(), len(cols),
                        'Check for %d columns in light curve FITS table' % len(
                            cols))
        self.test_value(table.nrows(), bins,
                        'Check for %d rows in light curve FITS table' % bins)
        for col in cols:
            self.test_assert(table.contains(col),
                             'FITS file contains "' + col + '" column')

        # Check EBOUNDS table
        table = fits['EBOUNDS']
        self._check_ebounds(table, nbins)

        # Close FITS file
        fits.close()

        # Return
        return

    def _check_output(self, filenameroot, bins, nreg, nout):
        """
        Check the output from a csphagen run
        """

        # OGIP files
        self._check_pha(filenameroot + '_pha_on.fits', bins)
        self._check_pha(filenameroot + '_pha_off.fits', bins)
        self._check_arf(filenameroot + '_arf.fits', bins)
        self._check_rmf(filenameroot + '_rmf.fits', bins)

        # Observations


        # Regions
        reg = gammalib.GSkyRegions(filenameroot + "_on.reg")
        self.test_value(reg.size(), 1, 'Check for ' + str(
            number) + ' region in source region file')
        reg = gammalib.GSkyRegions(filenameroot + "_off.reg")
        self.test_value(reg.size(), nreg, 'Check for ' + str(
            number) + ' region in background region file')

        return

    def _check_outobs(self, filenameroot, nout):
        """
        Check the output XML file containing ON/OFF observations
        """
        obs = gammalib.GObservations(filenameroot + '.xml')
        self.test_value(obs.size(), nout, 'Check for ' + str(
            number) + ' observations in XML file')

        return
