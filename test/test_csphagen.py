#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csphagen script.
#
# Copyright (C) 2017-2021 Luigi Tibaldo
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

        # Set test datasets and parameters
        self._myevents1      = self._datadir + '/crab_offaxis1.fits'
        self._myevents2      = self._datadir + '/crab_offaxis2.fits'
        self._exclusion      = self._datadir + '/crab_exclusion.fits'
        self._offevents      = self._datadir + '/background_events.fits'
        self._nreg_with_excl = 5
        self._nreg_wo_excl   = 8
        self._nreg_bkg_reg   = 5
        self._nreg_mul       = [self._nreg_with_excl, 6]
        self._nbins          = 10
        self._regfile_src    = self._datadir + '/crab_src_reg.reg'
        self._regfile_bkg    = self._datadir + '/crab_bkg_reg.reg'

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
        self.append(self._test_pickeling, 'Test csphagen pickeling')

        # Return
        return

    # Test csphagen on command line
    def _test_cmd(self):
        """
        Test csphagen on the command line.
        """
        # Set script name
        csphagen = self._script('csphagen')

        # Setup csphagen command
        cmd = csphagen + ' inobs="' + self._myevents1 + \
                         '" inmodel="NONE"' +\
                         ' caldb="' + self._caldb + '" irf="' + self._irf + \
                         '" ebinalg="LOG" emin=0.1 emax=100. enumbins=' + \
                         str(self._nbins) + \
                         ' etruemin=0.05 etruemax=150 etruebins=5'+ \
                         ' srcshape="CIRCLE" coordsys="CEL" ra=83.633 dec=22.0145' + \
                         ' rad=0.2 stack="no" inexclusion="' + \
                         self._exclusion + \
                         '" bkgmethod="REFLECTED" '+ \
                         'outobs="csphagen_cmd1_obs.xml" prefix="csphagen_cmd1" ' + \
                         'outmodel="csphagen_cmd1_model.xml" ' + \
                         'logfile="csphagen_cmd1.log" chatter=2'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
                         'Check successful execution from command line')

        # Check output files
        self._check_output('csphagen_cmd1', self._nbins, self._nreg_with_excl)
        self._check_outobs('csphagen_cmd1_obs.xml', 1)

        # Setup csphagen command
        cmd = csphagen + ' inobs="events_that_do_not_exist.fits" ' + \
                         'inmodel="NONE"' +\
                         ' caldb="' + self._caldb + '" irf="' + self._irf + \
                         '" ebinalg="LOG" emin=0.1 emax=100. enumbins=' + \
                         str(self._nbins) + \
                         ' etruemin=0.05 etruemax=150 etruebins=5'+ \
                         ' srcshape="CIRCLE" coordsys="CEL" ra=83.633 dec=22.0145' + \
                         ' rad=0.2 stack="no" inexclusion="' + \
                         self._exclusion + \
                         '" bkgmethod="REFLECTED" '+ \
                         'outobs="csphagen_cmd2_obs.xml" prefix="csphagen_cmd2" ' + \
                         'outmodel="csphagen_cmd2_model.xml" ' + \
                         'logfile="csphagen_cmd2.log" debug=yes chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check csphagen --help
        self._check_help(csphagen)

        # Return
        return

    def _test_python(self):
        """
        Test csphagen from Python
        """
        # Same test as from command line
        phagen = cscripts.csphagen()
        phagen['inobs']       = self._myevents1
        phagen['inmodel']     = 'NONE'
        phagen['caldb']       = self._caldb
        phagen['irf']         = self._irf
        phagen['ebinalg']     = 'LOG'
        phagen['emin']        = 0.1
        phagen['emax']        = 100.0
        phagen['enumbins']    = self._nbins
        phagen['srcshape']    = 'CIRCLE'
        phagen['coordsys']    = 'CEL'
        phagen['ra']          = 83.633
        phagen['dec']         = 22.0145
        phagen['rad']         = 0.2
        phagen['stack']       = False
        phagen['inexclusion'] = self._exclusion
        phagen['bkgmethod']   = 'REFLECTED'
        phagen['etruemin']    = 0.05
        phagen['etruemax']    = 150.0
        phagen['etruebins']   = 5
        phagen['outobs']      = 'csphagen_py1_obs.xml'
        phagen['outmodel']    = 'csphagen_py1_model.xml'
        phagen['prefix']      = 'csphagen_py1'
        phagen['logfile']     = 'csphagen_py1.log'
        phagen['chatter']     = 1

        # Execute script
        phagen.execute()

        # Check output
        self._check_output('csphagen_py1', self._nbins, self._nreg_with_excl)
        self._check_outobs('csphagen_py1_obs.xml', 1)

        # Now test without exclusion region
        phagen = cscripts.csphagen()
        phagen['inobs']     = self._myevents1
        phagen['inmodel']   = 'NONE'
        phagen['caldb']     = self._caldb
        phagen['irf']       = self._irf
        phagen['ebinalg']   = 'LOG'
        phagen['emin']      = 0.1
        phagen['emax']      = 100.0
        phagen['enumbins']  = self._nbins
        phagen['srcshape']  = 'CIRCLE'
        phagen['coordsys']  = 'CEL'
        phagen['ra']        = 83.633
        phagen['dec']       = 22.0145
        phagen['rad']       = 0.2
        phagen['stack']     = False
        phagen['bkgmethod'] = 'REFLECTED'
        phagen['etruemin']  = 0.05
        phagen['etruemax']  = 150.0
        phagen['etruebins'] = 5
        phagen['outobs']    = 'csphagen_py2_obs.xml'
        phagen['outmodel']  = 'csphagen_py2_model.xml'
        phagen['prefix']    = 'csphagen_py2'
        phagen['logfile']   = 'csphagen_py2.log'
        phagen['chatter']   = 2

        # Execute script
        phagen.execute()

        # Check output
        self._check_output('csphagen_py2', self._nbins, self._nreg_wo_excl)
        self._check_outobs('csphagen_py2_obs.xml', 1)

        # Test with multiple input observations, no stacking
        # Create observation container
        obs = gammalib.GObservations()
        for s, events in enumerate([self._myevents1, self._myevents2]):
            run = gammalib.GCTAObservation(events)
            run.id(str(s + 1))
            run.response(self._irf, gammalib.GCaldb('cta', self._caldb))
            obs.append(run)

        # Setup csphagen
        phagen = cscripts.csphagen(obs)
        phagen['inmodel']     = 'NONE'
        phagen['ebinalg']     = 'POW'
        phagen['emin']        = 0.1
        phagen['emax']        = 100.0
        phagen['enumbins']    = self._nbins
        phagen['ebingamma']   = 1.0
        phagen['srcshape']    = 'CIRCLE'
        phagen['coordsys']    = 'CEL'
        phagen['ra']          = 83.633
        phagen['dec']         = 22.0145
        phagen['rad']         = 0.2
        phagen['stack']       = False
        phagen['inexclusion'] = self._exclusion
        phagen['bkgmethod']   = 'REFLECTED'
        phagen['etruemin']    = 0.05
        phagen['etruemax']    = 150.0
        phagen['etruebins']   = 5
        phagen['outobs']      = 'csphagen_py3_obs.xml'
        phagen['outmodel']    = 'csphagen_py3_model.xml'
        phagen['prefix']      = 'csphagen_py3'
        phagen['logfile']     = 'csphagen_py3.log'
        phagen['chatter']     = 3

        # Run script
        phagen.execute()

        # Check output
        for s in range(2):
            self._check_output('csphagen_py3_' + str(s + 1), self._nbins,
                               self._nreg_mul[s])
        self._check_outobs('csphagen_py3_obs.xml', 2)

        # Setup csphagen for test with multiple input observations and stacking
        phagen = cscripts.csphagen(obs)
        phagen['inmodel']     = 'NONE'
        phagen['ebinalg']     = 'LOG'
        phagen['emin']        = 0.1
        phagen['emax']        = 100.0
        phagen['enumbins']    = self._nbins
        phagen['srcshape']    = 'CIRCLE'
        phagen['coordsys']    = 'CEL'
        phagen['ra']          = 83.633
        phagen['dec']         = 22.0145
        phagen['rad']         = 0.2
        phagen['stack']       = True
        phagen['inexclusion'] = self._exclusion
        phagen['bkgmethod']   = 'REFLECTED'
        phagen['etruemin']    = 0.05
        phagen['etruemax']    = 150.0
        phagen['etruebins']   = 5
        phagen['outobs']      = 'csphagen_py4_obs.xml'
        phagen['outmodel']    = 'csphagen_py4_model.xml'
        phagen['prefix']      = 'csphagen_py4'
        phagen['logfile']     = 'csphagen_py4.log'
        phagen['chatter']     = 4

        # Execute script
        phagen.execute()

        # Check output
        for s in range(2):
            self._check_output('csphagen_py4_stacked', self._nbins,
                               0, check_regions=False)
        self._check_outobs('csphagen_py4_obs.xml', 1)

        # Setup csphagen for test with custom On and Off regions provided
        phagen = cscripts.csphagen()
        phagen['inobs']      = self._myevents1
        phagen['inmodel']    = 'NONE'
        phagen['caldb']      = self._caldb
        phagen['irf']        = self._irf
        phagen['ebinalg']    = 'LOG'
        phagen['emin']       = 0.1
        phagen['emax']       = 100.0
        phagen['enumbins']   = self._nbins
        phagen['bkgmethod']  = 'CUSTOM'
        phagen['srcregfile'] = self._regfile_src
        phagen['bkgregfile'] = self._regfile_bkg
        phagen['etruemin']   = 0.05
        phagen['etruemax']   = 150.0
        phagen['etruebins']  = 5
        phagen['stack']      = False
        phagen['outobs']     = 'csphagen_py5_obs.xml'
        phagen['outmodel']   = 'csphagen_py5_model.xml'
        phagen['prefix']     = 'csphagen_py5'
        phagen['logfile']    = 'csphagen_py5.log'
        phagen['chatter']    = 2

        # Execute script
        phagen.execute()

        # Check output
        self._check_output('csphagen_py5', self._nbins, self._nreg_bkg_reg)
        self._check_outobs('csphagen_py5_obs.xml', 1)

        # Append off regions to observation container
        for run in obs:
            run.off_regions(gammalib.GSkyRegions(self._regfile_bkg))

        # Setup csphagen for test with multiple input observations and stacking
        phagen = cscripts.csphagen(obs)
        phagen['inmodel']     = 'NONE'
        phagen['ebinalg']     = 'LOG'
        phagen['emin']        = 0.1
        phagen['emax']        = 100.0
        phagen['enumbins']    = self._nbins
        phagen['coordsys']    = 'CEL'
        phagen['stack']       = True
        phagen['inexclusion'] = self._exclusion
        phagen['bkgmethod']   = 'CUSTOM'
        phagen['srcregfile']  = self._regfile_src
        phagen['etruemin']    = 0.05
        phagen['etruemax']    = 150.0
        phagen['etruebins']   = 5
        phagen['outobs']      = 'csphagen_py6_obs.xml'
        phagen['outmodel']    = 'csphagen_py6_model.xml'
        phagen['prefix']      = 'csphagen_py6'
        phagen['logfile']     = 'csphagen_py6.log'
        phagen['chatter']     = 4

        # Execute script
        phagen.execute()

        # Check output
        for s in range(2):
            self._check_output('csphagen_py6_stacked', self._nbins,
                               0, check_regions=False)
        self._check_outobs('csphagen_py6_obs.xml', 1)

        # Same test as from command line, but set exclusion with csphagen.exclusion_map
        phagen = cscripts.csphagen()
        phagen['inobs']       = self._myevents1
        phagen['inmodel']     = 'NONE'
        phagen['caldb']       = self._caldb
        phagen['irf']         = self._irf
        phagen['ebinalg']     = 'LOG'
        phagen['emin']        = 0.1
        phagen['emax']        = 100.0
        phagen['enumbins']    = self._nbins
        phagen['srcshape']    = 'CIRCLE'
        phagen['coordsys']    = 'CEL'
        phagen['ra']          = 83.633
        phagen['dec']         = 22.0145
        phagen['rad']         = 0.2
        phagen['stack']       = False
        phagen['bkgmethod']   = 'REFLECTED'
        phagen['etruemin']    = 0.05
        phagen['etruemax']    = 150.0
        phagen['etruebins']   = 5
        phagen['outobs']      = 'csphagen_py7_obs.xml'
        phagen['outmodel']    = 'csphagen_py7_model.xml'
        phagen['prefix']      = 'csphagen_py7'
        phagen['logfile']     = 'csphagen_py7.log'
        phagen['chatter']     = 1

        # Set exclusion map
        exclmap = gammalib.GSkyRegionMap(self._exclusion)
        phagen.exclusion_map(exclmap)

        # Execute script
        phagen.execute()

        # Check output
        self._check_output('csphagen_py7', self._nbins, self._nreg_with_excl)
        self._check_outobs('csphagen_py7_obs.xml', 1)

        # Test with rectangular region
        phagen = cscripts.csphagen()
        phagen['inobs']       = self._myevents1
        phagen['inmodel']     = 'NONE'
        phagen['caldb']       = self._caldb
        phagen['irf']         = self._irf
        phagen['ebinalg']     = 'LOG'
        phagen['emin']        = 0.1
        phagen['emax']        = 100.0
        phagen['enumbins']    = self._nbins
        phagen['srcshape']    = 'RECT'
        phagen['width']       = 0.15
        phagen['height']      = 0.25
        phagen['posang']      = 45.0
        phagen['coordsys']    = 'CEL'
        phagen['ra']          = 83.633
        phagen['dec']         = 22.0145
        phagen['rad']         = 0.2
        phagen['stack']       = False
        phagen['inexclusion'] = self._exclusion
        phagen['bkgmethod']   = 'REFLECTED'
        phagen['etruemin']    = 0.05
        phagen['etruemax']    = 150.0
        phagen['etruebins']   = 5
        phagen['outobs']      = 'csphagen_py8_obs.xml'
        phagen['outmodel']    = 'csphagen_py8_model.xml'
        phagen['prefix']      = 'csphagen_py8'
        phagen['logfile']     = 'csphagen_py8.log'
        phagen['chatter']     = 2

        # Execute script
        phagen.execute()

        # Check output
        self._check_output('csphagen_py8', self._nbins, 9)
        self._check_outobs('csphagen_py8_obs.xml', 1)

        # Setup csphagen for test with dedicated Off observations
        phagen = cscripts.csphagen()
        phagen['inobs']      = self._myevents1
        phagen['inobsoff']   = self._offevents
        phagen['inmodel']    = 'NONE'
        phagen['caldb']      = self._caldb
        phagen['irf']        = self._irf
        phagen['srcshape']   = 'CIRCLE'
        phagen['coordsys']   = 'CEL'
        phagen['ra']         = 83.633
        phagen['dec']        = 22.0145
        phagen['rad']        = 0.2
        phagen['ebinalg']    = 'LOG'
        phagen['emin']       = 0.1
        phagen['emax']       = 100.0
        phagen['enumbins']   = self._nbins
        phagen['bkgmethod']  = 'OFF'
        phagen['etruemin']   = 0.05
        phagen['etruemax']   = 150.0
        phagen['etruebins']  = 5
        phagen['stack']      = False
        phagen['outobs']     = 'csphagen_py9_obs.xml'
        phagen['outmodel']   = 'csphagen_py9_model.xml'
        phagen['prefix']     = 'csphagen_py9'
        phagen['logfile']    = 'csphagen_py9.log'
        phagen['chatter']    = 2

        # Execute script
        phagen.execute()

        # Check output
        self._check_output('csphagen_py9', self._nbins, 1)
        self._check_outobs('csphagen_py9_obs.xml', 1)

        # Setup csphagen for test with dedicated Off observations and
        # exclusion region
        phagen = cscripts.csphagen()
        phagen['inobs']       = self._myevents1
        phagen['inobsoff']    = self._offevents
        phagen['inmodel']     = 'NONE'
        phagen['caldb']       = self._caldb
        phagen['irf']         = self._irf
        phagen['srcshape']    = 'CIRCLE'
        phagen['coordsys']    = 'CEL'
        phagen['ra']          = 83.633
        phagen['dec']         = 22.0145
        phagen['rad']         = 0.2
        phagen['ebinalg']     = 'LOG'
        phagen['emin']        = 0.1
        phagen['emax']        = 100.0
        phagen['enumbins']    = self._nbins
        phagen['inexclusion'] = self._exclusion
        phagen['bkgmethod']   = 'OFF'
        phagen['etruemin']    = 0.05
        phagen['etruemax']    = 150.0
        phagen['etruebins']   = 5
        phagen['stack']       = False
        phagen['outobs']      = 'csphagen_py10_obs.xml'
        phagen['outmodel']    = 'csphagen_py10_model.xml'
        phagen['prefix']      = 'csphagen_py10'
        phagen['logfile']     = 'csphagen_py10.log'
        phagen['chatter']     = 2

        # Execute script
        phagen.execute()

        # Check output
        self._check_output('csphagen_py10', self._nbins, 1)
        self._check_outobs('csphagen_py10_obs.xml', 1)

        # Test with multiple input observations and dedicated Off runs
        # Create Off observation container
        obs_off = gammalib.GObservations()
        for s, events in enumerate([self._offevents, self._offevents]):
            run = gammalib.GCTAObservation(events)
            run.id(str(s + 1))
            run.response(self._irf, gammalib.GCaldb('cta', self._caldb))
            obs_off.append(run)

        # Setup csphagen
        phagen = cscripts.csphagen(obs)
        phagen.obs_off(obs_off)
        phagen['inmodel']    = 'NONE'
        phagen['srcshape']   = 'CIRCLE'
        phagen['coordsys']   = 'CEL'
        phagen['ra']         = 83.633
        phagen['dec']        = 22.0145
        phagen['rad']        = 0.2
        phagen['ebinalg']    = 'LOG'
        phagen['emin']       = 0.1
        phagen['emax']       = 100.0
        phagen['enumbins']   = self._nbins
        phagen['bkgmethod']  = 'OFF'
        phagen['etruemin']   = 0.05
        phagen['etruemax']   = 150.0
        phagen['etruebins']  = 5
        phagen['stack']      = False
        phagen['outobs']     = 'csphagen_py11_obs.xml'
        phagen['outmodel']   = 'csphagen_py11_model.xml'
        phagen['prefix']     = 'csphagen_py11'
        phagen['logfile']    = 'csphagen_py11.log'
        phagen['chatter']    = 2

        # Execute script
        phagen.execute()

        # Check output
        self._check_output('csphagen_py11_1', self._nbins, 1)
        self._check_output('csphagen_py11_2', self._nbins, 1)
        self._check_outobs('csphagen_py11_obs.xml', 2)

        # Return
        return

    # Test csphagen pickeling
    def _test_pickeling(self):
        """
        Test csphagen pickeling
        """
        # Perform pickeling tests of empty class
        self._pickeling(cscripts.csphagen())

        # Set-up csphagen
        phagen = cscripts.csphagen()
        phagen['inobs']       = self._myevents1
        phagen['inmodel']     = 'NONE'
        phagen['caldb']       = self._caldb
        phagen['irf']         = self._irf
        phagen['ebinalg']     = 'LOG'
        phagen['emin']        = 0.1
        phagen['emax']        = 100.0
        phagen['enumbins']    = self._nbins
        phagen['srcshape']    = 'CIRCLE'
        phagen['coordsys']    = 'CEL'
        phagen['ra']          = 83.633
        phagen['dec']         = 22.0145
        phagen['rad']         = 0.2
        phagen['stack']       = False
        phagen['inexclusion'] = self._exclusion
        phagen['bkgmethod']   = 'REFLECTED'
        phagen['etruemin']    = 0.05
        phagen['etruemax']    = 150.0
        phagen['etruebins']   = 5
        phagen['outobs']      = 'csphagen_py1_obs_pickle.xml'
        phagen['outmodel']    = 'csphagen_py1_model_pickle.xml'
        phagen['prefix']      = 'csphagen_py1_pickle'
        phagen['logfile']     = 'csphagen_py1_pickle.log'
        phagen['chatter']     = 1

        # Perform pickeling tests of filled class
        obj = self._pickeling(phagen)

        # Run csspec script and save light curve
        obj.logFileOpen()   # Make sure we get a log file
        obj.run()
        obj.save()

        # Check result files
        self._check_output('csphagen_py1_pickle', self._nbins, self._nreg_with_excl)
        self._check_outobs('csphagen_py1_obs_pickle.xml', 1)

        # Return
        return

    def _check_ebounds(self, table, bins):
        """
        Check EBOUNDS table
        """
        # Expected column names
        cols = ['E_MIN', 'E_MAX']

        # Check number of columns, rows and column names
        self.test_value(table.ncols(), len(cols), 'Check for %d columns in '
                        'energy bounds table' % len(cols))
        self.test_value(table.nrows(), bins,
                        'Check for %d rows in energy bounds table' % bins)
        for col in cols:
            self.test_assert(table.contains(col),
                             'Test if EBOUNDS contains "' + col + '" column')

        # Return
        return

    def _check_pha(self, filename, bins, ncols=8):
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
                        'Check for 3 extensions in PHA file')
        self.test_assert(fits.contains('SPECTRUM'),
                         'Test if PHA file contains "SPECTRUM" extension')
        self.test_assert(fits.contains('EBOUNDS'),
                         'Test if PHA file contains "EBOUNDS" extension')

        # Get SPECTRUM table
        table = fits['SPECTRUM']

        # Check FITS table structure
        self.test_value(table.ncols(), ncols,
                        'Check for %d columns in PHA table' % len(cols))
        self.test_value(table.nrows(), bins,
                        'Check for %d rows in PHA table' % bins)
        for col in cols:
            self.test_assert(table.contains(col),
                             'Test if PHA file contains "' + col + '" column')

        # Check EBOUNDS table
        table = fits['EBOUNDS']
        self._check_ebounds(table, bins)

        # Close FITS file
        fits.close()

        # Return
        return

    def _check_arf(self, filename, bins):
        """
        Check ARF file
        """
        # Expected column names
        cols = ['ENERG_LO', 'ENERG_HI', 'SPECRESP']

        # Open FITS file
        fits = gammalib.GFits(filename)

        # Check FITS file structure
        self.test_value(fits.size(), 2,
                        'Check for 2 extensions in ARF file')
        self.test_assert(fits.contains('SPECRESP'),
                         'Test if ARF file contains "SPECRESP" extension')

        # Get SPECRESP table
        table = fits['SPECRESP']

        # Check FITS table structure
        self.test_value(table.ncols(), len(cols),
                        'Check for %d columns in ARF table' % len(cols))
        self.test_value(table.nrows(), bins,
                        'Check for %d rows in ARF table' % bins)
        for col in cols:
            self.test_assert(table.contains(col),
                             'Test if ARF file contains "' + col + '" column')

        # Close FITS file
        fits.close()

        # Return
        return

    def _check_rmf(self, filename, bins, etruebins=17):
        """
        Check RMF file
        """
        # Expected column names
        cols = ['ENERG_LO', 'ENERG_HI', 'N_GRP', 'F_CHAN', 'N_CHAN', 'MATRIX']

        # Open FITS file
        fits = gammalib.GFits(filename)

        # Check FITS file structure
        self.test_value(fits.size(), 3,
                        'Check for 3 extensions in RMF file')
        self.test_assert(fits.contains('MATRIX'),
                         'Test if RMF file contains "MATRIX" extension')
        self.test_assert(fits.contains('EBOUNDS'),
                         'Test if RMF file contains "EBOUNDS" extension')

        # Get MATRIX table
        table = fits['MATRIX']

        # Check FITS table structure
        self.test_value(table.ncols(), len(cols),
                        'Check for %d columns in RMF table' % len(cols))
        self.test_value(table.nrows(), etruebins,
                        'Check for %d rows in RMF table' % etruebins)
        for col in cols:
            self.test_assert(table.contains(col),
                             'Test if RMF file contains "' + col + '" column')

        # Check EBOUNDS table
        table = fits['EBOUNDS']
        self._check_ebounds(table, bins)

        # Close FITS file
        fits.close()

        # Return
        return

    def _check_output(self, filenameroot, bins, nreg, check_regions=True):
        """
        Check the output from a csphagen run
        """
        # OGIP files
        self._check_pha(filenameroot + '_pha_on.fits', bins)
        self._check_pha(filenameroot + '_pha_off.fits', bins, ncols=9)
        self._check_arf(filenameroot + '_arf.fits', 17)
        self._check_rmf(filenameroot + '_rmf.fits', bins)

        # Optionally check for regions
        if check_regions:

            # Check On region
            onregion = (filenameroot+'_on.reg').replace('_1_', '_').replace('_2_', '_')
            reg      = gammalib.GSkyRegions(onregion)
            self.test_value(reg.size(), 1, 'Check for 1 region in source region file')

            # Check for Off region
            offregion = filenameroot + '_off.reg'
            reg       = gammalib.GSkyRegions(offregion)
            self.test_value(reg.size(), nreg, 'Check for ' + str(nreg) +
                            ' region in background region file')

        # Return
        return

    def _check_outobs(self, filenameroot, nout):
        """
        Check the output XML file containing ON/OFF observations
        """
        # Load observation container
        obs = gammalib.GObservations(filenameroot)

        # Check container size
        self.test_value(obs.size(), nout, 'Check for ' + str(nout) +
                        ' observations in XML file')

        # Return
        return

    def _read_pha_counts(self, filename):
        """
        Read and integrate the counts in a pha file.
        Pha file structure already tested in _check_pha()
        """
        # Open FITS file
        fits = gammalib.GFits(filename)

        # Get SPECTRUM table
        table = fits['SPECTRUM']

        # Integrate counts
        counts_col = table['COUNTS']
        counts     = 0
        for channel in range( counts_col.nrows() ):
            counts += counts_col[channel]

        # Close FITS file
        fits.close()

        # Return
        return counts
