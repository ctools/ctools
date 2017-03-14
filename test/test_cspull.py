#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cspull script.
#
# Copyright (C) 2016-2017 Juergen Knoedlseder
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
# Test class for cspull script #
# ============================ #
class Test(test):
    """
    Test class for cspull script

    This test class makes unit tests for the cspull script by using it
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
        self._inobs         = self._datadir + '/obs_stacked.xml'
        self._inobs_two     = self._datadir + '/obs_stacked_two.xml'
        self._stacked_model = self._datadir + '/crab_bkgcube.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('cspull')

        # Append tests
        self.append(self._test_cmd, 'Test cspull on command line')
        self.append(self._test_python, 'Test cspull from Python')

        # Return
        return

    # Test cspull on command line
    def _test_cmd(self):
        """
        Test cspull on the command line
        """
        # Set script name
        cspull = self._script('cspull')

        # Setup cspull command
        cmd = cspull+' inmodel="'+self._model+'"'+ \
                     ' outfile="cspull_cmd1.dat"'+ \
                     ' ntrials=2'+ \
                     ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                     ' ra=83.6331 dec=22.0145 emin=0.1 emax=100.0'+ \
                     ' enumbins=0 tmax=100.0 deadc=0.98 rad=5.0'+ \
                     ' npix=200 binsz=0.05'+ \
                     ' logfile="cspull_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check pull distribution file
        self._check_pull_file('cspull_cmd1.dat')

        # Setup cspull command
        cmd = cspull+' inmodel="model_that_does_not_exist.xml"'+ \
                     ' outfile="cspull_cmd1.dat"'+ \
                     ' ntrials=2'+ \
                     ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                     ' ra=83.6331 dec=22.0145 emin=0.1 emax=100.0'+ \
                     ' enumbins=0 tmax=100.0 deadc=0.98 rad=5.0'+ \
                     ' npix=200 binsz=0.05'+ \
                     ' logfile="cspull_cmd2.log"'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test cspull from Python
    def _test_python(self):
        """
        Test cspull from Python
        """
        # Set-up unbinned cspull
        pull = cscripts.cspull()
        pull['inmodel']  = self._model
        pull['outfile']  = 'cspull_py1.dat'
        pull['ntrials']  = 2
        pull['caldb']    = self._caldb
        pull['irf']      = self._irf
        pull['ra']       = 83.6331
        pull['dec']      = 22.0145
        pull['emin']     = 0.1
        pull['emax']     = 100.0
        pull['enumbins'] = 0
        pull['tmax']     = 100.0
        pull['deadc']    = 0.98
        pull['rad']      = 5.0
        pull['logfile']  = 'cspull_py1.log'
        pull['chatter']  = 2

        # Run cspull script
        pull.logFileOpen()   # Make sure we get a log file
        pull.run()
        #pull.save()

        # Check pull distribution file
        self._check_pull_file('cspull_py1.dat')

        # Set-up binned cspull
        pull = cscripts.cspull()
        pull['inmodel']  = self._model
        pull['outfile']  = 'cspull_py2.dat'
        pull['ntrials']  = 1
        pull['caldb']    = self._caldb
        pull['irf']      = self._irf
        pull['ra']       = 83.6331
        pull['dec']      = 22.0145
        pull['emin']     = 0.1
        pull['emax']     = 100.0
        pull['enumbins'] = 10
        pull['tmax']     = 100.0
        pull['deadc']    = 0.98
        pull['rad']      = 5.0
        pull['npix']     = 20
        pull['binsz']    = 0.2
        pull['coordsys'] = 'CEL'
        pull['proj']     = 'TAN'
        pull['logfile']  = 'cspull_py2.log'
        pull['chatter']  = 3

        # Execute cspull script
        pull.execute()

        # Check pull distribution file
        self._check_pull_file('cspull_py2.dat', rows=2)

        # Set-up cspull from event list
        pull = cscripts.cspull()
        pull['inobs']    = self._events
        pull['inmodel']  = self._model
        pull['outfile']  = 'cspull_py3.dat'
        pull['ntrials']  = 1
        pull['caldb']    = self._caldb
        pull['irf']      = self._irf
        pull['enumbins'] = 0
        pull['logfile']  = 'cspull_py3.log'
        pull['chatter']  = 4

        # Execute cspull script
        pull.execute()

        # Check pull distribution file
        self._check_pull_file('cspull_py3.dat', rows=2)

        # Build observation container with unbinned observation
        cta = gammalib.GCTAObservation(self._events)
        obs = gammalib.GObservations()
        obs.append(cta)

        # Set-up cspull from observation container with unbinned observation
        pull = cscripts.cspull(obs)
        pull['inmodel']  = self._model
        pull['outfile']  = 'cspull_py4.dat'
        pull['ntrials']  = 1
        pull['caldb']    = self._caldb
        pull['irf']      = self._irf
        pull['enumbins'] = 0
        pull['logfile']  = 'cspull_py4.log'
        pull['chatter']  = 4

        # Execute cspull script
        pull.execute()

        # Check pull distribution file
        self._check_pull_file('cspull_py4.dat', rows=2)

        # Set-up stacked cspull with one observation where response cubes
        # are specified in the input observation
        pull = cscripts.cspull()
        pull['inobs']    = self._inobs
        pull['inmodel']  = self._stacked_model
        pull['outfile']  = 'cspull_py5.dat'
        pull['ntrials']  = 1
        pull['enumbins'] = 0
        pull['logfile']  = 'cspull_py5.log'
        pull['chatter']  = 4

        # Execute cspull script
        pull.execute()

        # Check pull distribution file
        self._check_pull_file('cspull_py5.dat', rows=2)

        # Set-up stacked cspull with two observations from which response
        # cubes will be computed internally. The IRFs are also specified
        # in the input observation, hence we do not need to query the IRF
        # parameters.
        pull = cscripts.cspull()
        pull['inobs']    = self._inobs_two
        pull['inmodel']  = self._model
        pull['outfile']  = 'cspull_py6.dat'
        pull['ntrials']  = 1
        pull['emin']     = 0.02
        pull['emax']     = 100.0
        pull['enumbins'] = 10
        pull['npix']     = 20
        pull['binsz']    = 0.2
        pull['coordsys'] = 'CEL'
        pull['proj']     = 'TAN'
        pull['logfile']  = 'cspull_py6.log'
        pull['chatter']  = 4

        # Execute cspull script
        pull.execute()

        # Check pull distribution file
        self._check_pull_file('cspull_py6.dat', rows=2)

        # Return
        return

    # Check pull file
    def _check_pull_file(self, filename, rows=3):
        """
        Check pull file
        """
        # Open pull file as CSV file
        pulls = gammalib.GCsv(filename)

        # Check dimensions
        self.test_value(pulls.nrows(), rows,
                        'Check for number of rows in pull file')

        # Return
        return
