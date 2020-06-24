#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csscs script.
#
# Copyright (C) 2020 Luigi Tibaldo
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
import os
import gammalib
import cscripts
from testing import test


# =========================== #
# Test class for csscs script #
# =========================== #
class Test(test):
    """
    Test class for csscs script

    This test class makes unit tests for the csscs script by using it from
    the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set test datasets and parameters
        self._evs_offaxis = self._datadir + '/crab_offaxis1.fits'
        self._exclusion   = self._datadir + '/crab_exclusion.fits'
        self._expcube     = self._datadir + '/crab_expcube.fits'
        self._psfcube     = self._datadir + '/crab_psfcube.fits'
        self._bkgcube     = self._datadir + '/crab_bkgcube.fits'
        self._bkgmodel    = self._datadir + '/crab_bkgcube.xml'
        self._edispcube   = self._datadir + '/crab_edispcube.fits'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csscs')

        # Append tests
        self.append(self._test_cmd, 'Test csscs on command line')
        self.append(self._test_python, 'Test csscs from Python')
        self.append(self._test_pickeling, 'Test csscs pickeling')

        # Return
        return

    # Test csscs on command line
    def _test_cmd(self):
        """
        Test csscs on the command line
        """
        # Set script name
        csscs = self._script('csscs')

        # Setup csscs command
        cmd = csscs + ' inobs="' + self._cntcube + '"' + \
              ' expcube="' + self._expcube + '"' + \
              ' psfcube="' + self._psfcube + '"' + \
              ' bkgcube="' + self._bkgcube + '"' + \
              ' inmodel="' + self._bkgmodel + '"' + \
              ' srcnames="Crab"' + \
              ' emin=1. emax=100.0 nxpix=2 nypix=2' + \
              ' binsz=0.1 coordsys="CEL" proj="CAR"' + \
              ' xref=83.63 yref=22.01 rad=0.2' + \
              ' outfile="csscs_cmd1.fits"' + \
              ' logfile="csscs_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_value(self._execute(cmd), 0,
                        'Check successful execution from command line')

        # Setup csscs command
        cmd = csscs + ' inobs="' + self._cntcube + '"' + \
              ' expcube="' + self._expcube + '"' + \
              ' psfcube="' + self._psfcube + '"' + \
              ' bkgcube="' + self._bkgcube + '"' + \
              ' inmodel="' + self._bkgmodel + '"' + \
              ' srcnames="Crab;source_that_is_not_in_the_model"' + \
              ' emin=1. emax=100.0 nxpix=2 nypix=2' + \
              ' binsz=0.1 coordsys="CEL" proj="CAR"' + \
              ' xref=83.63 yref=22.01 rad=0.2' + \
              ' outfile="csscs_cmd1.fits"' + \
              ' logfile="csscs_cmd1.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
                         'Check invalid source name when executed from command line')

        # Check csscs --help
        self._check_help(csscs)

        # Return
        return

    # Test csscs from Python
    def _test_python(self):
        """
        Test csscs from Python
        """
        # Allocate csscs
        script = cscripts.csscs()

        # Check that saving does not nothing
        script['logfile'] = 'csscs_py0.log'
        script['outfile'] = 'csscs_py0.fits'
        script.logFileOpen()
        script.save()
        self.test_assert(not os.path.isfile('csscs_py0.fits'),
                         'Check that no FITS file has been created')

        # Check that clearing does not lead to an exception or segfault
        script.clear()

        # Set up csscs for simple binned analysis
        script['logfile'] = 'csscs_py1.log'
        script['inobs'] = self._cntcube
        script['expcube'] = self._expcube
        script['psfcube'] = self._psfcube
        script['bkgcube'] = self._bkgcube
        script['inmodel'] = self._bkgmodel
        script['srcnames'] = 'Crab'
        script['emin'] = 1.
        script['emax'] = 100.
        script['nxpix'] = 2
        script['nypix'] = 2
        script['binsz'] = 0.1
        script['coordsys'] = 'CEL'
        script['proj'] = 'CAR'
        script['xref'] = 83.63
        script['yref'] = 22.01
        script['rad'] = 0.2
        script['outfile'] = 'csscs_py1.fits'
        script['chatter'] = 4

        # Run csscs script and save results
        script.logFileOpen()
        script.execute()

        # Check result
        self._check_result_file('csscs_py1.fits')

        # Set up for same test without TS computation
        script = cscripts.csscs()
        script['logfile']  = 'csscs_py2.log'
        script['inobs']    = self._cntcube
        script['expcube']  = self._expcube
        script['psfcube']  = self._psfcube
        script['bkgcube']  = self._bkgcube
        script['inmodel']  = self._bkgmodel
        script['srcnames'] = 'Crab'
        script['emin']     = 1.
        script['emax']     = 100.
        script['nxpix']    = 2
        script['nypix']    = 2
        script['binsz']    = 0.1
        script['coordsys'] = 'CEL'
        script['proj']     = 'CAR'
        script['xref']     = 83.63
        script['yref']     = 22.01
        script['rad']      = 0.2
        script['calc_ts']  = False
        script['outfile']  = 'csscs_py2.fits'
        script['chatter']  = 2

        # Run csscs script and save results
        script.logFileOpen()
        script.execute()

        # Check result
        self._check_result_file('csscs_py2.fits', calc_ts=False)

        # # Set up for same test with upper limit computation
        script = cscripts.csscs()
        script['logfile']   = 'csscs_py3.log'
        script['inobs']     = self._cntcube
        script['expcube']   = self._expcube
        script['psfcube']   = self._psfcube
        script['bkgcube']   = self._bkgcube
        script['inmodel']   = self._bkgmodel
        script['srcnames']  = 'Crab'
        script['emin']      = 1.
        script['emax']      = 100.
        script['nxpix']     = 2
        script['nypix']     = 2
        script['binsz']     = 0.1
        script['coordsys']  = 'CEL'
        script['proj']      = 'CAR'
        script['xref']      = 83.63
        script['yref']      = 22.01
        script['rad']       = 0.2
        script['calc_ts']   = True
        script['calc_ulim'] = True
        script['outfile']   = 'csscs_py3.fits'
        script['chatter']   = 2

        # Run csscs script and save results
        script.logFileOpen()
        script.execute()

        # Check result
        self._check_result_file('csscs_py3.fits', calc_ulim=True)

        # # Check that results can be piped in memory
        self._check_results_mem(script, calc_ulim=True)

        # Set up for test without multiprocessing
        script = cscripts.csscs()
        script['logfile']  = 'csscs_py4.log'
        script['inobs']    = self._cntcube
        script['expcube']  = self._expcube
        script['psfcube']  = self._psfcube
        script['bkgcube']  = self._bkgcube
        script['inmodel']  = self._bkgmodel
        script['srcnames'] = 'Crab'
        script['emin']     = 1.
        script['emax']     = 100.
        script['nxpix']    = 2
        script['nypix']    = 2
        script['binsz']    = 0.1
        script['coordsys'] = 'CEL'
        script['proj']     = 'CAR'
        script['xref']     = 83.63
        script['yref']     = 22.01
        script['rad']      = 0.2
        script['outfile']  = 'csscs_py4.fits'
        script['nthreads'] = 1
        script['chatter']  = 2

        # Run csscs script and save results
        script.logFileOpen()
        script.execute()

        # Check result
        self._check_result_file('csscs_py4.fits')

        # Set up for test with energy dispersion
        script = cscripts.csscs()
        script['logfile']   = 'csscs_py5.log'
        script['inobs']     = self._cntcube
        script['expcube']   = self._expcube
        script['psfcube']   = self._psfcube
        script['bkgcube']   = self._bkgcube
        script['edisp']     = True
        script['edispcube'] = self._edispcube
        script['inmodel']   = self._bkgmodel
        script['srcnames']  = 'Crab'
        script['emin']      = 1.
        script['emax']      = 100.
        script['nxpix']     = 2
        script['nypix']     = 2
        script['binsz']     = 0.1
        script['coordsys']  = 'CEL'
        script['proj']      = 'CAR'
        script['xref']      = 83.63
        script['yref']      = 22.01
        script['rad']       = 0.2
        script['outfile']   = 'csscs_py5.fits'
        script['chatter']   = 2

        # Run csscs script and save results
        script.logFileOpen()
        script.execute()

        # Check result
        self._check_result_file('csscs_py5.fits')

        # Set up for test in On/Off mode with background model
        script = cscripts.csscs()
        script['logfile']     = 'csscs_py6.log'
        script['inobs']       = self._evs_offaxis
        script['caldb']       = self._caldb
        script['irf']         = self._irf
        script['inmodel']     = self._model
        script['srcnames']    = 'Crab'
        script['emin']        = 1.
        script['emax']        = 100.
        script['nxpix']       = 2
        script['nypix']       = 2
        script['binsz']       = 0.1
        script['coordsys']    = 'CEL'
        script['proj']        = 'CAR'
        script['xref']        = 83.63
        script['yref']        = 22.01
        script['rad']         = 0.2
        script['method']      = 'ONOFF'
        script['inexclusion'] = self._exclusion
        script['enumbins']    = 2
        script['outfile']     = 'csscs_py6.fits'
        script['chatter']     = 2

        # Run csscs script and save results
        script.logFileOpen()
        script.execute()

        # Check result
        self._check_result_file('csscs_py6.fits')

        # Same test with exclusion map piped in memory
        script = cscripts.csscs()
        script['logfile']  = 'csscs_py7.log'
        script['inobs']    = self._evs_offaxis
        script['caldb']    = self._caldb
        script['irf']      = self._irf
        script['inmodel']  = self._model
        script['srcnames'] = 'Crab'
        script['emin']     = 1.
        script['emax']     = 100.
        script['nxpix']    = 2
        script['nypix']    = 2
        script['binsz']    = 0.1
        script['coordsys'] = 'CEL'
        script['proj']     = 'CAR'
        script['xref']     = 83.63
        script['yref']     = 22.01
        script['rad']      = 0.2
        script['method']   = 'ONOFF'
        script['enumbins'] = 2
        script['outfile']  = 'csscs_py7.fits'
        script['chatter']  = 2

        # Set exclusion map
        exclmap = gammalib.GSkyRegionMap(self._exclusion)
        script.exclusion_map(exclmap)

        # Run csscs script and save results
        script.logFileOpen()
        script.execute()

        # Check result
        self._check_result_file('csscs_py7.fits')

        # Set up for test in On/Off mode without background model
        script = cscripts.csscs()
        script['logfile']       = 'csscs_py8.log'
        script['inobs']         = self._evs_offaxis
        script['caldb']         = self._caldb
        script['irf']           = self._irf
        script['inmodel']       = self._model
        script['srcnames']      = 'Crab'
        script['emin']          = 1.
        script['emax']          = 100.
        script['nxpix']         = 2
        script['nypix']         = 2
        script['binsz']         = 0.1
        script['coordsys']      = 'CEL'
        script['proj']          = 'CAR'
        script['xref']          = 83.63
        script['yref']          = 22.01
        script['rad']           = 0.2
        script['method']        = 'ONOFF'
        script['use_model_bkg'] = False
        script['inexclusion']   = self._exclusion
        script['enumbins']      = 2
        script['outfile']       = 'csscs_py8.fits'
        script['chatter']       = 2

        # Run csscs script and save results
        script.logFileOpen()
        script.execute()

        # Check result
        self._check_result_file('csscs_py8.fits')

        # Set up for test in unbinned mode
        script = cscripts.csscs()
        script['logfile']  = 'csscs_py9.log'
        script['inobs']    = self._evs_offaxis
        script['caldb']    = self._caldb
        script['irf']      = self._irf
        script['inmodel']  = self._model
        script['srcnames'] = 'Crab'
        script['emin']     = 1.
        script['emax']     = 100.
        script['nxpix']    = 2
        script['nypix']    = 2
        script['binsz']    = 0.1
        script['coordsys'] = 'CEL'
        script['proj']     = 'CAR'
        script['xref']     = 83.63
        script['yref']     = 22.01
        script['rad']      = 0.2
        script['method']   = 'UNBINNED'
        script['enumbins'] = 2
        script['outfile']  = 'csscs_py9.fits'
        script['chatter']  = 2

        # Run csscs script and save results
        script.logFileOpen()
        script.execute()

        # Check result
        self._check_result_file('csscs_py9.fits')

        # Return
        return

    # Test csscs pickeling
    def _test_pickeling(self):
        """
        Test csscs pickeling
        """
        # Perform pickeling test of empty class
        self._pickeling(cscripts.csscs())

        # Setup script for pickling text
        script = cscripts.csscs()
        script['logfile']  = 'csscs_py1_pickle.log'
        script['inobs']    = self._cntcube
        script['expcube']  = self._expcube
        script['psfcube']  = self._psfcube
        script['bkgcube']  = self._bkgcube
        script['inmodel']  = self._bkgmodel
        script['srcnames'] = 'Crab'
        script['emin']     = 1.
        script['emax']     = 100.
        script['nxpix']    = 2
        script['nypix']    = 2
        script['binsz']    = 0.1
        script['coordsys'] = 'CEL'
        script['proj']     = 'CAR'
        script['xref']     = 83.63
        script['yref']     = 22.01
        script['rad']      = 0.2
        script['outfile']  = 'csscs_py1_pickle.fits'
        script['chatter']  = 2

        # Perform pickeling tests of filled class
        obj = self._pickeling(script)

        # Run csscs script and save result
        obj.logFileOpen()  # Make sure we get a log file
        obj.execute()

        # Check result file
        self._check_result_file('csscs_py1_pickle.fits')

        # Return
        return

    # Check access to results in memory'
    def _check_results_mem(self, script, srcnames=['Crab'], nx=2, ny=2, calc_ts=True,
                           calc_ulim=False):
        """
        Check results in memory

        Parameters
        ----------
        script : `~cscripts.csscs`
            csscs object
        srcnames : list, optional
            Source names
        nx : int, optional
            Number of X pixels
        ny : int, optional
            Number of Y pixels
        calc_ts : bool, optional
            Check TS HDU
        calc_ulim : bool, optional
            Check upper limit HDU
        """
        # Check that the entire fits container can be extracted
        fits = script.fits()

        # Dump to disk and verify
        fits.saveto('csscs_py3_mem.fits', True)
        self._check_result_file('csscs_py3_mem.fits', srcnames=srcnames, nx=nx, ny=ny,
                                calc_ts=calc_ts, calc_ulim=calc_ulim)

        # Now check direct access in memory for skymaps
        # Loop over sources
        for name in srcnames:
            # Check flux map
            skymap = script.flux(name)
            self._check_map(skymap, nx=nx, ny=ny)
            # Check flux error map
            skymap = script.flux_error(name)
            self._check_map(skymap, nx=nx, ny=ny)
            # If requested check TS map
            if calc_ts:
                skymap = script.ts(name)
                self._check_map(skymap, nx=nx, ny=ny)
            # If requested check upper limit map
            if calc_ulim:
                skymap = script.ulimit(name)
                self._check_map(skymap, nx=nx, ny=ny)

        return

    # Check result file
    def _check_result_file(self, filename, srcnames=['Crab'], nx=2, ny=2,
                           calc_ts=True, calc_ulim=False):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Sky map file name
        srcnames : list, optional
            Source names
        nx : int, optional
            Number of X pixels
        ny : int, optional
            Number of Y pixels
        calc_ts : bool, optional
            Check TS HDU
        calc_ulim : bool, optional
            Check upper limit HDU
        """
        # Loop over source names
        for name in srcnames:
            # Check flux map
            skymap = gammalib.GSkyMap(filename + '[' + name + ' FLUX]')
            self._check_map(skymap, nx=nx, ny=ny)
            # Check flux error map
            skymap = gammalib.GSkyMap(filename + '[' + name + ' FLUX ERROR]')
            self._check_map(skymap, nx=nx, ny=ny)
            # If requested check TS map
            if calc_ts:
                skymap = gammalib.GSkyMap(filename + '[' + name + ' TS]')
                self._check_map(skymap, nx=nx, ny=ny)
            # If requested check upper limit map
            if calc_ulim:
                skymap = gammalib.GSkyMap(filename + '[' + name + ' FLUX UPPER LIMIT]')
                self._check_map(skymap, nx=nx, ny=ny)

        # Return
        return

    # Check sky map
    def _check_map(self, skymap, nx=2, ny=2):
        """
        Check sky map

        Parameters
        ----------
        skymap : `~gammalib.GSkyMap`
            Sky map
        nx : int, optional
            Number of X pixels
        ny : int, optional
            Number of Y pixels
        """
        # Check dimensions
        self.test_value(skymap.nmaps(), 1, 'Check number of maps')
        self.test_value(skymap.nx(), nx, 'Check for number of X pixels')
        self.test_value(skymap.ny(), ny, 'Check for number of Y pixels')

        # Return
        return
