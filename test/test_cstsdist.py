#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the cstsdist script.
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


# ============================== #
# Test class for cstsdist script #
# ============================== #
class Test(test):
    """
    Test class for cstsdist script

    This test class makes unit tests for the cstsdist script by using it
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
        self._stacked_inobs = self._datadir + '/obs_stacked.xml'
        self._stacked_model = self._datadir + '/crab_bkgcube.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('cstsdist')

        # Append tests
        self.append(self._test_cmd, 'Test cstsdist on command line')
        self.append(self._test_python_unbinned,
                    'Test cstsdist from Python in unbinned mode')
        self.append(self._test_python_binned,
                    'Test cstsdist from Python in binned mode')
        self.append(self._test_python_stacked,
                    'Test cstsdist from Python in stacked mode')
        self.append(self._test_pickeling, 'Test cstsdist pickeling')

        # Return
        return

    # Test cstsdist on command line
    def _test_cmd(self):
        """
        Test cstsdist on the command line
        """
        # Set script name
        cstsdist = self._script('cstsdist')

        # Setup cstsdist command
        cmd = cstsdist+' inobs="'+self._events+'"'+ \
                       ' inmodel="'+self._model+'" srcname="Crab"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                       ' ntrials=2 '+ \
                       ' outfile="cstsdist_cmd1.dat"'+ \
                       ' logfile="cstsdist_cmd1.log" chatter=2'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('cstsdist_cmd1.dat')

        # Setup cstsdist command
        cmd = cstsdist+' inobs="event_file_that_does_not_exist.xml"'+ \
                       ' inmodel="'+self._model+'" srcname="Crab"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                       ' ntrials=1 '+ \
                       ' outfile="cstsdist_cmd2.dat"'+ \
                       ' logfile="cstsdist_cmd2.log" chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check cstsdist --help
        self._check_help(cstsdist)

        # Return
        return

    # Test cstsdist from Python
    def _test_python_unbinned(self):
        """
        Test cstsdist from Python
        """
        # Set-up cstsdist for event list
        tsdist = cscripts.cstsdist()
        tsdist['inobs']   = self._events
        tsdist['inmodel'] = self._model
        tsdist['srcname'] = 'Crab'
        tsdist['caldb']   = self._caldb
        tsdist['irf']     = self._irf
        tsdist['ntrials'] = 2
        tsdist['outfile'] = 'cstsdist_py1.dat'
        tsdist['logfile'] = 'cstsdist_py1.log'
        tsdist['chatter'] = 3

        # Run cstsdist script
        tsdist.logFileOpen()   # Make sure we get a log file
        tsdist.run()
        tsdist.save()

        # Check TS distribution file
        self._check_result_file('cstsdist_py1.dat')

        # Set-up cstsdist without multiprocessing
        tsdist = cscripts.cstsdist()
        tsdist['inobs']    = self._events
        tsdist['inmodel']  = self._model
        tsdist['srcname']  = 'Crab'
        tsdist['caldb']    = self._caldb
        tsdist['irf']      = self._irf
        tsdist['ntrials']  = 2
        tsdist['outfile']  = 'cstsdist_py2.dat'
        tsdist['logfile']  = 'cstsdist_py2.log'
        tsdist['chatter']  = 3
        tsdist['nthreads'] = 1

        # Run cstsdist script
        tsdist.logFileOpen()  # Make sure we get a log file
        tsdist.run()
        tsdist.save()

        # Check TS distribution file
        self._check_result_file('cstsdist_py2.dat')

        # Return
        return

    def _test_python_binned(self):
        """
        Test cstsdist from Python
        """
        # Set-up cstsdist for counts cube
        tsdist = cscripts.cstsdist()
        tsdist['inobs']   = self._cntcube
        tsdist['inmodel'] = self._model
        tsdist['srcname'] = 'Crab'
        tsdist['expcube'] = 'NONE'
        tsdist['psfcube'] = 'NONE'
        tsdist['bkgcube'] = 'NONE'
        tsdist['caldb']   = self._caldb
        tsdist['irf']     = self._irf
        tsdist['ntrials'] = 2
        tsdist['outfile'] = 'cstsdist_py2.dat'
        tsdist['logfile'] = 'cstsdist_py2.log'
        tsdist['chatter'] = 4

        # Run cstsdist script
        tsdist.logFileOpen()   # Make sure we get a log file
        tsdist.execute()

        # Check TS distribution file
        self._check_result_file('cstsdist_py2.dat')

        # Return
        return

    def _test_python_stacked(self):
        """
        Test cstsdist from Python
        """
        # Set-up cstsdist for stacked cube
        tsdist = cscripts.cstsdist()
        tsdist['inobs']   = self._stacked_inobs
        tsdist['inmodel'] = self._stacked_model
        tsdist['srcname'] = 'Crab'
        tsdist['ntrials'] = 2
        tsdist['outfile'] = 'cstsdist_py3.dat'
        tsdist['logfile'] = 'cstsdist_py3.log'
        tsdist['chatter'] = 4

        # Run cstsdist script
        tsdist.logFileOpen()   # Make sure we get a log file
        tsdist.execute()

        # Check TS distribution file
        self._check_result_file('cstsdist_py3.dat')

        # Return
        return

    # Test cstsdist pickeling
    def _test_pickeling(self):
        """
        Test cstsdist pickeling
        """
        # Perform pickeling tests of empty class
        self._pickeling(cscripts.cstsdist())

        # Set-up cstsdist for event list
        tsdist = cscripts.cstsdist()
        tsdist['inobs']   = self._events
        tsdist['inmodel'] = self._model
        tsdist['srcname'] = 'Crab'
        tsdist['caldb']   = self._caldb
        tsdist['irf']     = self._irf
        tsdist['ntrials'] = 2
        tsdist['outfile'] = 'cstsdist_py1_pickle.dat'
        tsdist['logfile'] = 'cstsdist_py1_pickle.log'
        tsdist['chatter'] = 3

        # Perform pickeling tests of filled class
        obj = self._pickeling(tsdist)

        # Run csphasecrv script and save light curve
        obj.logFileOpen()   # Make sure we get a log file
        obj.run()
        obj.save()

        # Check TS distribution file
        self._check_result_file('cstsdist_py1_pickle.dat')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nrows=3, ncols=9):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Name of result file
        nrows : int, optional
            Required number of rows
        ncols : int, optional
            Required number of columns
        """
        # Open result file as CSV file
        results = gammalib.GCsv(filename, ',')

        # Check dimensions
        self.test_value(results.nrows(), nrows, 'Check rows in TS file')
        self.test_value(results.ncols(), ncols, 'Check columns in TS file')

        # Return
        return
