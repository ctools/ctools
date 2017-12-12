#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csphasecrv script.
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
import gammalib
import cscripts
from testing import test


# ================================ #
# Test class for csphasecrv script #
# ================================ #
class Test(test):
    """
    Test class for csphasecrv script

    This test class makes unit tests for the csphasecrv script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set test data
        self._phased_events  = self._datadir + '/phased_events.fits'
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
        self.name('csphasecrv')

        # Append tests
        self.append(self._test_cmd, 'Test csphasecrv on command line')
        self.append(self._test_python, 'Test csphasecrv from Python')

        # Return
        return

    # Test csphasecrv on command line
    def _test_cmd(self):
        """
        Test csphasecrv on the command line.
        """
        # Set script name
        csphasecrv = self._script('csphasecrv')

        # Setup csphasecrv command
        cmd = csphasecrv+' inobs="'+self._phased_events+'"'+ \
                         ' inmodel="'+self._model+'" srcname="Crab"'+ \
                         ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                         ' phbinalg="LIN" phbins=5 method="3D"'+ \
                         ' enumbins=0 emin=0.1 emax=100.0'+ \
                         ' outfile="csphasecrv_cmd1.fits"'+ \
                         ' logfile="csphasecrv_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check phase curve
        self._check_phase_curve('csphasecrv_cmd1.fits', 5)

        # Setup csphasecrv command
        cmd = csphasecrv+' inobs="events_that_do_not_exist.fits"'+ \
                         ' inmodel="'+self._model+'" srcname="Crab"'+ \
                         ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                         ' phbinalg="LIN" phbins=5 method="3D"'+ \
                         ' enumbins=0 emin=0.1 emax=100.0'+ \
                         ' outfile="csphasecrv_cmd2.fits"'+ \
                         ' logfile="csphasecrv_cmd2.log" debug=yes chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Check csphasecrv --help
        self._check_help(csphasecrv)

        # Return
        return

    # Test csphasecrv from Python
    def _test_python(self):
        """
        Test csphasecrv from Python
        """
        # Set-up unbinned csphasecrv
        pcrv = cscripts.csphasecrv()
        pcrv['inobs']    = self._phased_events
        pcrv['inmodel']  = self._model
        pcrv['srcname']  = 'Crab'
        pcrv['caldb']    = self._caldb
        pcrv['irf']      = self._irf
        pcrv['phbinalg'] = 'LIN'
        pcrv['phbins']   = 5
        pcrv['method']   = '3D'
        pcrv['enumbins'] = 0
        pcrv['emin']     = 0.1
        pcrv['emax']     = 100.0
        pcrv['outfile']  = 'csphasecrv_py1.fits'
        pcrv['logfile']  = 'csphasecrv_py1.log'
        pcrv['chatter']  = 2
        pcrv['publish']  = True

        # Run csphasecrv script and save phase curve
        pcrv.logFileOpen()   # Make sure we get a log file
        pcrv.run()
        pcrv.save()

        # Check phase curve
        self._check_phase_curve('csphasecrv_py1.fits', 5)

        # Now use FILE as time bin algorithm. For this we need first to
        # create an ASCII file. We use now 2 phase bins. The ASCII file
        # is saved into the file "csphasecrv_py2.dat".
        csv     = gammalib.GCsv(2,2)
        phmin   = 0.0
        phdelta = 1.0 / float(csv.nrows())
        for i in range(csv.nrows()):
            csv[i,0] = '%.5f' % (phmin +  i    * phdelta)
            csv[i,1] = '%.5f' % (phmin + (i+1) * phdelta)
        csv.save('csphasecrv_py2.dat', ' ', True)

        # Set-up unbinned csphasecrv
        pcrv = cscripts.csphasecrv()
        pcrv['inobs']     = self._phased_events
        pcrv['inmodel']   = self._model
        pcrv['srcname']   = 'Crab'
        pcrv['caldb']     = self._caldb
        pcrv['irf']       = self._irf
        pcrv['phbinalg']  = 'FILE'
        pcrv['phbinfile'] = 'csphasecrv_py2.dat'
        pcrv['method']    = '3D'
        pcrv['enumbins']  = 0
        pcrv['emin']      = 0.1
        pcrv['emax']      = 100.0
        pcrv['outfile']   = 'csphasecrv_py2.fits'
        pcrv['logfile']   = 'csphasecrv_py2.log'
        pcrv['chatter']   = 3

        # Execute csphasecrv script
        pcrv.logFileOpen()   # Make sure we get a log file
        pcrv.execute()

        # Check phase curve
        self._check_phase_curve('csphasecrv_py2.fits', csv.nrows())

        # Now we setup an observation container on input. We attached the
        # model to the observation container so that csphasecrv should
        # no longer query for the parameter.
        cta = gammalib.GCTAObservation(self._phased_events)
        obs = gammalib.GObservations()
        obs.append(cta)
        obs.models(self._model)

        # Set-up unbinned csphasecrv from observation container. Now use
        # the GTI algorithm so that we test all timing algorithms.
        pcrv = cscripts.csphasecrv(obs)
        pcrv['srcname']  = 'Crab'
        pcrv['caldb']    = self._caldb
        pcrv['irf']      = self._irf
        pcrv['phbinalg'] = 'LIN'
        pcrv['phbins']   = 5
        pcrv['method']   = '3D'
        pcrv['enumbins'] = 0
        pcrv['emin']     = 0.1
        pcrv['emax']     = 100.0
        pcrv['outfile']  = 'csphasecrv_py3.fits'
        pcrv['logfile']  = 'csphasecrv_py3.log'
        pcrv['chatter']  = 4

        # Execute csphasecrv script
        pcrv.logFileOpen()   # Make sure we get a log file
        pcrv.execute()

        # Check phase curve
        self._check_phase_curve('csphasecrv_py3.fits', 5)

        # Binned csphasecrv
        pcrv = cscripts.csphasecrv()
        pcrv['inobs']    = self._phased_events
        pcrv['inmodel']  = self._model
        pcrv['srcname']  = 'Crab'
        pcrv['caldb']    = self._caldb
        pcrv['irf']      = self._irf
        pcrv['phbinalg'] = 'LIN'
        pcrv['phbins']   = 2
        pcrv['method']   = '3D'
        pcrv['emin']     = 0.1
        pcrv['emax']     = 100.0
        pcrv['enumbins'] = 10
        pcrv['coordsys'] = 'CEL'
        pcrv['proj']     = 'TAN'
        pcrv['xref']     = 83.63
        pcrv['yref']     = 22.01
        pcrv['nxpix']    = 20
        pcrv['nypix']    = 20
        pcrv['binsz']    = 0.02
        pcrv['outfile']  = 'csphasecrv_py4.fits'
        pcrv['logfile']  = 'csphasecrv_py4.log'
        pcrv['chatter']  = 4

        # Execute csphasecrv script
        pcrv.logFileOpen()   # Make sure we get a log file
        pcrv.execute()

        # Check phase curve
        self._check_phase_curve('csphasecrv_py4.fits', 2)

        # On/Off csphasecrv
        pcrv = cscripts.csphasecrv()
        pcrv['inobs']     = self._offaxis_events
        pcrv['inmodel']   = self._model_onoff
        pcrv['srcname']   = 'Crab'
        pcrv['caldb']     = self._caldb
        pcrv['irf']       = self._irf
        pcrv['phbinalg']  = 'LIN'
        pcrv['phbins']    = 2
        pcrv['method']    = 'ONOFF'
        pcrv['statistic'] = 'WSTAT'
        pcrv['emin']      = 0.1
        pcrv['emax']      = 100.0
        pcrv['enumbins']  = 10
        pcrv['coordsys']  = 'CEL'
        pcrv['xref']      = 83.63
        pcrv['yref']      = 22.01
        pcrv['rad']       = 0.2
        pcrv['etruemin']  = 0.05
        pcrv['etruemax']  = 150.0
        pcrv['etruebins'] = 5
        pcrv['outfile']   = 'csphasecrv_py5.fits'
        pcrv['logfile']   = 'csphasecrv_py5.log'
        pcrv['chatter']   = 4

        # Execute csphasecrv script
        pcrv.logFileOpen()   # Make sure we get a log file
        pcrv.execute()

        # Check phase curve
        self._check_phase_curve('csphasecrv_py5.fits', 2)

        # Return
        return

    # Check phase curve result file
    def _check_phase_curve(self, filename, bins):
        """
        Check phase curve file
        """
        # Expected column names
        cols = ['PHASE_MIN', 'PHASE_MAX', 'Prefactor', 'e_Prefactor',
                'Index', 'e_Index']

        # Open FITS file
        fits = gammalib.GFits(filename)

        # Check FITS file structure
        self.test_value(fits.size(), 2,
             'Check for 2 extensions in phase curve FITS file')
        self.test_assert(fits.contains('PHASECURVE'),
             'FITS file contains "PHASECURVE" extension')

        # Get PHASECURVE table
        table = fits['PHASECURVE']

        # Check FITS table structure
        self.test_value(table.ncols(), len(cols),
             'Check for %d columns in phase curve FITS table' % len(cols))
        self.test_value(table.nrows(), bins,
             'Check for %d rows in phase curve FITS table' % bins)
        for col in cols:
            self.test_assert(table.contains(col),
                 'FITS file contains "'+col+'" column')

        # Check that table has been filled
        # Prefactor has right order of magnitude
        for s in range(table.nrows()):
            self.test_value(table['Prefactor'][s], 4.e-16, 3.e-16,
                            'Check prefactor value')

        # Close FITS file
        fits.close()
        
        # Return
        return
