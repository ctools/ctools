#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csspec script.
#
# Copyright (C) 2016-2018 Juergen Knoedlseder
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


# ============================ #
# Test class for csspec script #
# ============================ #
class Test(test):
    """
    Test class for csspec script

    This test class makes unit tests for the csspec script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set members
        self._inonoff     = self._datadir + '/onoff_obs.xml'
        self._expcube     = self._datadir + '/crab_expcube.fits'
        self._psfcube     = self._datadir + '/crab_psfcube.fits'
        self._edispcube   = self._datadir + '/crab_edispcube.fits'
        self._bkgcube     = self._datadir + '/crab_bkgcube.fits'
        self._bkgmodel    = self._datadir + '/crab_bkgcube.xml'
        self._onoff_model = self._datadir + '/onoff_model.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csspec')

        # Append tests
        self.append(self._test_cmd, 'Test csspec on command line')
        self.append(self._test_python, 'Test csspec from Python')
        self.append(self._test_pickeling, 'Test csspec pickeling')

        # Return
        return

    # Test csspec on command line
    def _test_cmd(self):
        """
        Test csspec on the command line
        """
        # Set script name
        csspec = self._script('csspec')

        # Setup csspec command
        cmd = csspec+' inobs="'+self._events+'"'+ \
                     ' inmodel="'+self._model+'"'+ \
                     ' srcname="Crab" method="AUTO"'+ \
                     ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                     ' ebinalg="LOG" enumbins=2 emin=1.0 emax=100.0'+ \
                     ' outfile="csspec_cmd1.fits"'+ \
                     ' logfile="csspec_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('csspec_cmd1.fits', 2)

        # Setup csspec command
        cmd = csspec+' inobs="input_file_that_does_not_exist.fits"'+ \
                     ' inmodel="'+self._model+'"'+ \
                     ' srcname="Crab" method="AUTO"'+ \
                     ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                     ' ebinalg="LOG" enumbins=2 emin=1.0 emax=100.0'+ \
                     ' outfile="csspec_cmd2.fits"'+ \
                     ' logfile="csspec_cmd2.log" debug=yes chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check csspec --help
        self._check_help(csspec)

        # Return
        return

    # Test csspec from Python
    def _test_python(self):
        """
        Test csspec from Python
        """
        # Set-up unbinned csspec (no multi-processing)
        spec = cscripts.csspec()
        spec['inobs']    = self._events
        spec['inmodel']  = self._model
        spec['srcname']  = 'Crab'
        spec['caldb']    = self._caldb
        spec['irf']      = self._irf
        spec['method']   = 'AUTO'
        spec['ebinalg']  = 'LOG'
        spec['enumbins'] = 2
        spec['emin']     = 1.0
        spec['emax']     = 100.0
        spec['nthreads'] = 1
        spec['outfile']  = 'csspec_py1.fits'
        spec['logfile']  = 'csspec_py1.log'
        spec['chatter']  = 2
        spec['publish']  = True

        # Run csspec script
        spec.logFileOpen()   # Make sure we get a log file
        spec.run()
        spec.save()

        # Check result file
        self._check_result_file('csspec_py1.fits', 2)

        # Set-up stacked csspec (no multi-processing)
        spec = cscripts.csspec()
        spec['inobs']     = self._cntcube
        spec['expcube']   = self._expcube
        spec['psfcube']   = self._psfcube
        spec['edispcube'] = self._edispcube
        spec['bkgcube']   = self._bkgcube
        spec['inmodel']   = self._bkgmodel
        spec['srcname']   = 'Crab'
        spec['method']    = 'AUTO'
        spec['ebinalg']   = 'LOG'
        spec['enumbins']  = 2
        spec['emin']      = 0.1
        spec['emax']      = 100.0
        spec['nthreads']  = 1
        spec['outfile']   = 'csspec_py2.fits'
        spec['logfile']   = 'csspec_py2.log'
        spec['chatter']   = 3
        spec['publish']   = False

        # Execute csspec script
        spec.execute()

        # Check result file
        self._check_result_file('csspec_py2.fits', 2)

        # Set-up On/Off csspec (no multi-processing)
        spec = cscripts.csspec()
        spec['inobs']    = self._inonoff
        spec['inmodel']  = self._onoff_model
        spec['srcname']  = 'Crab'
        spec['method']   = 'AUTO'
        spec['ebinalg']  = 'LOG'
        spec['enumbins'] = 2
        spec['emin']     = 1.0
        spec['emax']     = 10.0
        spec['nthreads'] = 1
        spec['outfile']  = 'csspec_py3.fits'
        spec['logfile']  = 'csspec_py3.log'
        spec['chatter']  = 4
        spec['publish']  = False

        # Execute csspec script
        spec.execute()

        # Check result file
        self._check_result_file('csspec_py3.fits', 2)

        # Set-up On/Off csspec with NODES method
        spec = cscripts.csspec()
        spec['inobs']    = self._inonoff
        spec['inmodel']  = self._onoff_model
        spec['srcname']  = 'Crab'
        spec['method']   = 'NODES'
        spec['ebinalg']  = 'LOG'
        spec['enumbins'] = 2
        spec['emin']     = 1.0
        spec['emax']     = 10.0
        spec['outfile']  = 'csspec_py4.fits'
        spec['logfile']  = 'csspec_py4.log'
        spec['chatter']  = 4
        spec['publish']  = False

        # Execute csspec script
        spec.execute()

        # Check result file
        self._check_result_file('csspec_py4.fits', 2)

        # Set-up On/Off csspec using WSTAT
        spec = cscripts.csspec()
        spec['inobs']     = self._inonoff
        spec['inmodel']   = self._onoff_model
        spec['srcname']   = 'Crab'
        spec['method']    = 'AUTO'
        spec['statistic'] = 'WSTAT'
        spec['ebinalg']   = 'LOG'
        spec['enumbins']  = 2
        spec['emin']      = 1.0
        spec['emax']      = 10.0
        spec['outfile']   = 'csspec_py5.fits'
        spec['logfile']   = 'csspec_py5.log'
        spec['chatter']   = 2

        # Execute csspec script
        spec.execute()

        # Check result file
        self._check_result_file('csspec_py5.fits', 2)

        # Set-up unbinned csspec with likelihood profile (no multi-processing)
        spec = cscripts.csspec()
        spec['inobs']       = self._events
        spec['inmodel']     = self._model
        spec['srcname']     = 'Crab'
        spec['caldb']       = self._caldb
        spec['irf']         = self._irf
        spec['method']      = 'AUTO'
        spec['ebinalg']     = 'LOG'
        spec['enumbins']    = 2
        spec['emin']        = 1.0
        spec['emax']        = 100.0
        spec['nthreads']    = 1
        spec['dll_sigstep'] = 2.5
        spec['dll_sigmax']  = 5
        spec['outfile']     = 'csspec_py6.fits'
        spec['logfile']     = 'csspec_py6.log'
        spec['chatter']     = 2
        spec['publish']     = True

        # Run csspec script
        spec.execute()

        # Check result file
        self._check_result_file('csspec_py6.fits', 2, 5)

        # Return
        return

    # Test csspec pickeling
    def _test_pickeling(self):
        """
        Test csspec pickeling
        """
        # Perform pickeling tests of empty class
        self._pickeling(cscripts.csspec())

        # Set-up unbinned csspec
        spec = cscripts.csspec()
        spec['inobs']    = self._events
        spec['inmodel']  = self._model
        spec['srcname']  = 'Crab'
        spec['caldb']    = self._caldb
        spec['irf']      = self._irf
        spec['method']   = 'AUTO'
        spec['ebinalg']  = 'LOG'
        spec['enumbins'] = 2
        spec['emin']     = 1.0
        spec['emax']     = 100.0
        spec['outfile']  = 'csspec_py1_pickle.fits'
        spec['logfile']  = 'csspec_py1_pickle.log'
        spec['chatter']  = 2
        spec['publish']  = True

        # Perform pickeling tests of filled class
        obj = self._pickeling(spec)

        # Run csspec script and save light curve
        obj.logFileOpen()   # Make sure we get a log file
        obj.run()
        obj.save()

        # Check result file
        self._check_result_file('csspec_py1_pickle.fits', 2)

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, bins, lbins=0):
        """
        Check result file
        """
        # Open result file
        fits = gammalib.GFits(filename)

        # Get spectrum table
        spectrum = fits['SPECTRUM']

        # Check dimensions
        self.test_value(spectrum.nrows(), bins,
             'Check for %d rows in spectrum' % bins)

        # Make sure that the header keywords are present
        req_keywords = ['SED_TYPE', 'UL_CONF']
        for keyword in req_keywords:
            self.test_assert(spectrum.has_card(keyword), 
                'Check that SED has \'%s\' keyword' % keyword)

        # Do checks specific to likelihood profile table
        flux_keys = spectrum.card('SED_TYPE').string()
        flux_keys = flux_keys.split(',')
        if lbins > 0:
            # Should be 13 columns
            self.test_value(spectrum.ncols(), 13,
                'Check for 8 columns in spectrum')
            self.test_value(flux_keys[0], 'likelihood',
                'Check that \'SED_TYPE\' card is \'likelihood\'')
            # Make sure the profile has the appropriate entries
            self.test_value(spectrum['norm_scan'].elements(0), lbins,
                'Check for %d rows in flux scan profile' % lbins)
            self.test_value(spectrum['dloglike_scan'].elements(0), lbins,
                'Check for %d rows in likelihood profile' % lbins)
        else:
            # Should be 10 columns
            self.test_value(spectrum.ncols(), 10,
                'Check for 8 columns in spectrum')
            self.test_value(flux_keys[0], 'norm',
                'Check that \'SED_TYPE\' card is \'norm\'')
        
        # Return
        return
