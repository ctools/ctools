#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cslightcrv script.
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


# ================================ #
# Test class for cslightcrv script #
# ================================ #
class Test(test):
    """
    Test class for cslightcrv script

    This test class makes unit tests for the cslightcrv script by using it
    from the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Add off-axis events for classical analysis
        self._offaxis_events = self._datadir + '/crab_offaxis1.fits'
        self._model_onoff    = self._datadir + '/crab_onoff.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('cslightcrv')

        # Append tests
        self.append(self._test_cmd, 'Test cslightcrv on command line')
        self.append(self._test_python, 'Test cslightcrv from Python')
        self.append(self._test_pickeling, 'Test cslightcrv pickeling')

        # Return
        return

    # Test cslightcrv on command line
    def _test_cmd(self):
        """
        Test cslightcrv on the command line.
        """
        # Set script name
        cslightcrv = self._script('cslightcrv')

        # Setup cslightcrv command
        cmd = cslightcrv+' inobs="'+self._events+'"'+ \
                         ' inmodel="'+self._model+'" srcname="Crab"'+ \
                         ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                         ' tbinalg=LIN '+ \
                         ' tmin=2020-01-01T00:00:00 tmax=2020-01-01T00:05:00'+ \
                         ' tbins=3 method=3D enumbins=0 emin=0.1 emax=100.0'+ \
                         ' outfile="lightcurve_cmd1.fits"'+ \
                         ' logfile="cslightcrv_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check light curve
        self._check_light_curve('lightcurve_cmd1.fits', 3)

        # Setup cslightcrv command
        cmd = cslightcrv+' inobs="events_that_do_not_exist.fits"'+ \
                         ' inmodel="'+self._model+'" srcname="Crab"'+ \
                         ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                         ' tbinalg=LIN '+ \
                         ' tmin=2020-01-01T00:00:00 tmax=2020-01-01T00:05:00'+ \
                         ' tbins=3 method=3D enumbins=0 emin=0.1 emax=100.0'+ \
                         ' outfile="lightcurve_cmd1.fits"'+ \
                         ' logfile="cslightcrv_cmd2.log" debug=yes'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check cslightcrv --help
        self._check_help(cslightcrv)

        # Return
        return

    # Test cslightcrv from Python
    def _test_python(self):
        """
        Test cslightcrv from Python
        """
        # Set-up unbinned cslightcrv
        lcrv = cscripts.cslightcrv()
        lcrv['inobs']    = self._events
        lcrv['inmodel']  = self._model
        lcrv['srcname']  = 'Crab'
        lcrv['caldb']    = self._caldb
        lcrv['irf']      = self._irf
        lcrv['tbinalg']  = 'LIN'
        lcrv['tmin']     = '2020-01-01T00:00:00'
        lcrv['tmax']     = '2020-01-01T00:05:00'
        lcrv['tbins']    = 3
        lcrv['method']   = '3D'
        lcrv['enumbins'] = 0
        lcrv['emin']     = 0.1
        lcrv['emax']     = 100.0
        lcrv['outfile']  = 'cslightcrv_py1.fits'
        lcrv['logfile']  = 'cslightcrv_py1.log'
        lcrv['chatter']  = 2
        lcrv['publish']  = True

        # Run cslightcrv script and save light curve
        lcrv.logFileOpen()   # Make sure we get a log file
        lcrv.run()
        lcrv.save()

        # Check light curve
        self._check_light_curve('cslightcrv_py1.fits', 3)

        # Now use FILE as time bin algorithm. For this we need first to
        # create an ASCII file. We use now 6 time bins. The ASCII file
        # is saved into the file "lightcurve_py2.dat".
        csv    = gammalib.GCsv(2,2)
        tmin   = 58849.00
        tdelta = 0.0017361
        for i in range(csv.nrows()):
            csv[i,0] = '%.5f' % (tmin +  i    * tdelta)
            csv[i,1] = '%.5f' % (tmin + (i+1) * tdelta)
        csv.save('cslightcrv_py2.dat', ' ', True)

        # Set-up unbinned cslightcrv
        lcrv = cscripts.cslightcrv()
        lcrv['inobs']    = self._events
        lcrv['inmodel']  = self._model
        lcrv['srcname']  = 'Crab'
        lcrv['caldb']    = self._caldb
        lcrv['irf']      = self._irf
        lcrv['tbinalg']  = 'FILE'
        lcrv['tbinfile'] = 'cslightcrv_py2.dat'
        lcrv['method']   = '3D'
        lcrv['enumbins'] = 0
        lcrv['emin']     = 0.1
        lcrv['emax']     = 100.0
        lcrv['fix_bkg']  =  True
        lcrv['outfile']  = 'cslightcrv_py2.fits'
        lcrv['logfile']  = 'cslightcrv_py2.log'
        lcrv['chatter']  = 3

        # Execute cslightcrv script
        lcrv.execute()

        # Check light curve
        self._check_light_curve('cslightcrv_py2.fits', 2)

        # Now we setup an observation container on input. We attached the
        # model to the observation container so that cslightcrv should
        # no longer query for the parameter.
        cta = gammalib.GCTAObservation(self._events)
        obs = gammalib.GObservations()
        obs.append(cta)
        obs.models(self._model)

        # Set-up unbinned cslightcrv from observation container. Now use
        # the GTI algorithm so that we test all timing algorithms.
        lcrv = cscripts.cslightcrv(obs)
        lcrv['srcname']  = 'Crab'
        lcrv['caldb']    = self._caldb
        lcrv['irf']      = self._irf
        lcrv['tbinalg']  = 'GTI'
        lcrv['method']   = '3D'
        lcrv['enumbins'] = 0
        lcrv['emin']     = 0.1
        lcrv['emax']     = 100.0
        lcrv['outfile']  = 'cslightcrv_py3.fits'
        lcrv['logfile']  = 'cslightcrv_py3.log'
        lcrv['chatter']  = 4

        # Execute cslightcrv script
        lcrv.execute()

        # Check light curve
        self._check_light_curve('cslightcrv_py3.fits', 1)

        # Binned cslightcrv
        lcrv = cscripts.cslightcrv()
        lcrv['inobs']    = self._events
        lcrv['inmodel']  = self._model
        lcrv['srcname']  = 'Crab'
        lcrv['caldb']    = self._caldb
        lcrv['irf']      = self._irf
        lcrv['tbinalg']  = 'LIN'
        lcrv['tmin']     = '2020-01-01T00:00:00'
        lcrv['tmax']     = '2020-01-01T00:05:00'
        lcrv['tbins']    = 2
        lcrv['method']   = '3D'
        lcrv['emin']     = 0.1
        lcrv['emax']     = 100.0
        lcrv['enumbins'] = 10
        lcrv['coordsys'] = 'CEL'
        lcrv['proj']     = 'TAN'
        lcrv['xref']     = 83.63
        lcrv['yref']     = 22.01
        lcrv['nxpix']    = 20
        lcrv['nypix']    = 20
        lcrv['binsz']    = 0.02
        lcrv['outfile']  = 'cslightcrv_py4.fits'
        lcrv['logfile']  = 'cslightcrv_py4.log'
        lcrv['chatter']  = 4

        # Execute cslightcrv script
        lcrv.execute()

        # Check light curve
        self._check_light_curve('cslightcrv_py4.fits', 2)

        # cslightcrv with classical analysis
        lcrv = cscripts.cslightcrv()
        lcrv['inobs']     = self._offaxis_events
        lcrv['inmodel']   = self._model_onoff
        lcrv['srcname']   = 'Crab'
        lcrv['caldb']     = self._caldb
        lcrv['irf']       = self._irf
        lcrv['tbinalg']   = 'LIN'
        lcrv['tmin']      = '2020-01-01T00:00:00'
        lcrv['tmax']      = '2020-01-01T00:05:00'
        lcrv['tbins']     = 2
        lcrv['method']    = 'ONOFF'
        lcrv['emin']      = 0.1
        lcrv['emax']      = 100.0
        lcrv['enumbins']  = 10
        lcrv['coordsys']  = 'CEL'
        lcrv['xref']      = 83.63
        lcrv['yref']      = 22.01
        lcrv['rad']       = 0.2
        lcrv['etruemin']  = 0.05
        lcrv['etruemax']  = 150.0
        lcrv['etruebins'] = 20
        lcrv['statistic'] = 'WSTAT'
        lcrv['outfile']   = 'cslightcrv_py5.fits'
        lcrv['logfile']   = 'cslightcrv_py5.log'
        lcrv['chatter']   = 4

        # Execute cslightcrv script
        lcrv.execute()

        # Check light curve
        self._check_light_curve('cslightcrv_py5.fits', 2)

        # Set-up cslightcrv without multiprocessing
        lcrv = cscripts.cslightcrv()
        lcrv['inobs'] = self._events
        lcrv['inmodel'] = self._model
        lcrv['srcname'] = 'Crab'
        lcrv['caldb'] = self._caldb
        lcrv['irf'] = self._irf
        lcrv['tbinalg'] = 'LIN'
        lcrv['tmin'] = '2020-01-01T00:00:00'
        lcrv['tmax'] = '2020-01-01T00:05:00'
        lcrv['tbins'] = 3
        lcrv['method'] = '3D'
        lcrv['enumbins'] = 0
        lcrv['emin'] = 0.1
        lcrv['emax'] = 100.0
        lcrv['outfile'] = 'cslightcrv_py6.fits'
        lcrv['logfile'] = 'cslightcrv_py6.log'
        lcrv['chatter'] = 2
        lcrv['publish'] = True
        lcrv['nthreads'] = 1

        # Run cslightcrv script and save light curve
        lcrv.logFileOpen()  # Make sure we get a log file
        lcrv.run()
        lcrv.save()

        # Check light curve
        self._check_light_curve('cslightcrv_py6.fits', 3)

        # Return
        return

    # Test cslightcrv pickeling
    def _test_pickeling(self):
        """
        Test cslightcrv pickeling
        """
        # Perform pickeling tests of empty class
        self._pickeling(cscripts.cslightcrv())

        # Set-up unbinned cslightcrv
        lcrv = cscripts.cslightcrv()
        lcrv['inobs']    = self._events
        lcrv['inmodel']  = self._model
        lcrv['srcname']  = 'Crab'
        lcrv['caldb']    = self._caldb
        lcrv['irf']      = self._irf
        lcrv['tbinalg']  = 'LIN'
        lcrv['tmin']     = '2020-01-01T00:00:00'
        lcrv['tmax']     = '2020-01-01T00:05:00'
        lcrv['tbins']    = 3
        lcrv['method']   = '3D'
        lcrv['enumbins'] = 0
        lcrv['emin']     = 0.1
        lcrv['emax']     = 100.0
        lcrv['outfile']  = 'cslightcrv_py1_pickle.fits'
        lcrv['logfile']  = 'cslightcrv_py1_pickle.log'
        lcrv['chatter']  = 2
        lcrv['publish']  = True

        # Perform pickeling tests of filled class
        obj = self._pickeling(lcrv)

        # Run cslightcrv script and save light curve
        obj.logFileOpen()   # Make sure we get a log file
        obj.run()
        obj.save()

        # Check light curve
        self._check_light_curve('cslightcrv_py1_pickle.fits', 3)

        # Return
        return

    # Check light curve result file
    def _check_light_curve(self, filename, bins, prefactor=5.7e-16):
        """
        Check light curve file
        """
        # Expected column names
        cols = ['MJD', 'e_MJD', 'Prefactor', 'e_Prefactor',
                'Index', 'e_Index', 'TS', 'DiffUpperLimit',
                'FluxUpperLimit', 'EFluxUpperLimit']

        # Open FITS file
        fits = gammalib.GFits(filename)

        # Check FITS file structure
        self.test_value(fits.size(), 2,
             'Check for 2 extensions in light curve FITS file')
        self.test_assert(fits.contains('LIGHTCURVE'),
             'FITS file contains "LIGHTCURVE" extension')

        # Get LIGHTCURVE table
        table = fits['LIGHTCURVE']

        # Check FITS table structure
        self.test_value(table.ncols(), len(cols),
             'Check for %d columns in light curve FITS table' % len(cols))
        self.test_value(table.nrows(), bins,
             'Check for %d rows in light curve FITS table' % bins)
        for col in cols:
            self.test_assert(table.contains(col),
                 'FITS file contains "'+col+'" column')

        # Check that the Prefactor has the right order of magnitude
        for s in range(table.nrows()):
            self.test_value(table['Prefactor'][s], prefactor, 0.2*prefactor,
                            'Check prefactor value')

        # Close FITS file
        fits.close()

        # Return
        return
