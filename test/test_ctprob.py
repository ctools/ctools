#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctprob tool.
#
# Copyright (C) 2017 by Leonardo Di Venere
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
import ctools
from testing import test


# ============================ #
# Test class for ctprob tool #
# ============================ #
class Test(test):
    """
    Test class for ctprob tool
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set test data
        self._invalid_events = self._datadir + '/invalid_event_list.fits'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('ctprob')

        # Append tests
        self.append(self._test_cmd, 'Test ctprob on command line')
        self.append(self._test_python, 'Test ctprob from Python')

        # Return
        return

    # Test ctprob on command line
    def _test_cmd(self):
        """
        Test ctprob on the command line
        """
        # Set tool name
        ctprob = self._tool('ctprob')

        # Setup ctprob command
        cmd = ctprob+' inobs="'+self._events+'"'+ \
                     ' outobs="ctprob_cmd1.fits"'+ \
                     ' inmodel='+self._model+ \
                     ' caldb='+self._caldb+' irf='+self._irf+ \
                     ' edisp=no'+ \
                     ' logfile="ctprob_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctprob_cmd1.fits')

        # Setup ctprob command
        cmd = ctprob+' inobs="event_file_that_does_not_exist.fits"'+ \
                     ' outobs="ctprob_cmd2.fits"'+ \
                     ' inmodel='+self._model+ \
                     ' caldb='+self._caldb+' irf='+self._irf+ \
                     ' edisp=no'+ \
                     ' logfile="ctprob_cmd2.log" debug=yes chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Check ctprob --help
        self._check_help(ctprob)

        # Return
        return

    # Test ctprob from Python
    def _test_python(self):
        """
        Test ctprob from Python
        """
        # Allocate empty ctprob tool
        prob = ctools.ctprob()

        # Check that empty ctprob tool holds an empty observation
        self._check_obs(prob.obs(), nobs=0)

        # Check that saving does nothing
        prob['outobs']  = 'ctprob_py0.fits'
        prob['logfile'] = 'ctprob_py0.log'
        prob.logFileOpen()
        prob.save()
        self.test_assert(not os.path.isfile('ctprob_py0.fits'),
             'Check that no event list has been created')

        # Check that clearing does not lead to an exception or segfault
        prob.clear()

        # Now set ctprob parameters
        prob['inobs']   = self._events
        prob['outobs']  = 'ctprob_py1.fits'
        prob['inmodel'] = self._model
        prob['caldb']   = self._caldb
        prob['irf']     = self._irf
        prob['edisp']   = False
        prob['publish'] = True
        prob['logfile'] = 'ctprob_py1.log'
        prob['chatter'] = 2

        # Run ctprob tool
        prob.logFileOpen()   # Make sure we get a log file
        prob.run()
        prob.save()

        # Check result file
        self._check_result_file('ctprob_py1.fits')

        # Setup ctprob tool for an invalid event file
        prob = ctools.ctprob()
        prob['inobs']   = self._invalid_events
        prob['outobs']  = 'ctprob_py2.fits'
        prob['inmodel'] = self._model
        prob['caldb']   = self._caldb
        prob['irf']     = self._irf
        prob['edisp']   = False
        prob['logfile'] = 'ctprob_py2.log'
        prob['chatter'] = 3

        # Execute tool
        prob.logFileOpen()  # Needed to get a new log file
        self.test_try('Test invalid event file')
        try:
            prob.execute()
            self.test_try_failure('Exception not thrown')
        except RuntimeError:
            self.test_try_success()

        # Setup ctprob with a list of observations and test prefix
        prob = ctools.ctprob()
        prob['inobs']   = self._datadir + '/obs_unbinned.xml'
        prob['prefix']  = 'ctprob_py3_'
        prob['outobs']  = 'ctprob_py3_obs.xml'
        prob['inmodel'] = self._model
        prob['caldb']   = self._caldb
        prob['irf']     = self._irf
        prob['edisp']   = True
        prob['logfile'] = 'ctprob_py3.log'
        prob['chatter'] = 4

        # Run ctprob tool
        prob.logFileOpen()   # Make sure we get a log file
        prob.execute()

        # Check result file
        self._check_result_file('ctprob_py3_crab_events.fits')

        # Copy ctprob tool
        cpy_prob = prob.copy()

        # Clear phase tool
        prob.clear()

        # Execute copy
        cpy_prob['prefix']  = 'ctprob_py4_'
        cpy_prob['outobs']  = 'ctprob_py4_obs.xml'
        cpy_prob['logfile'] = 'ctprob_py4.log'
        cpy_prob['chatter'] = 4
        cpy_prob.logFileOpen()  # Needed to get a new log file
        cpy_prob.execute()

        # Check result file
        self._check_result_file('ctprob_py4_ctprob_py3_crab_events.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename, nevents=3736):
        """
        Check result file

        Parameters
        ----------
        filename : str
            Event list file name
        nevents : int, optional
            Expected number of events
        """
        # Open result file
        events = gammalib.GCTAEventList(filename)

        # Check event list
        self._check_events(events, nevents=nevents)

        # Check that probability columns exist
        colname1 = 'PROB_CTABackgroundModel'
        colname2 = 'PROB_Crab'
        self._check_column(filename, colname1)
        self._check_column(filename, colname2)        

        # Check that the probability values are correctly normalized
        self._check_normalization(filename, [colname1,colname2])

        # Return
        return

    # Check observation and event list
    def _check_obs(self, obs, nobs=1, nevents=22220):
        """
        Check observation and event list

        Parameters
        ----------
        obs : `~gammalib.GObservations`
            Models
        nobs : int, optional
            Expected number of observations
        nevents : int, optional
            Expected number of events
        """
        # Check number of observations
        self.test_value(obs.size(), nobs, 'Check number of observations')

        # If there is an observation then check the event list of the first
        # one
        if obs.size() > 0:
            self._check_events(obs[0].events(), nevents=nevents)

        # Return
        return

    # Check events
    def _check_events(self, events, nevents=22220):
        """
        Check event list

        Parameters
        ----------
        events : `~gammalib.GCTAEventList`
            Event list
        nevents : int, optional
            Expected number of events
        """
        # Check event list
        self.test_value(events.size(), nevents, 'Check number of events')

        # Return
        return


    def _check_column(self, filename, colname):
        """
        Check that column colname exists.

        Parameters
        ----------
        filename : str
            Event list file name
        colname : str
            Name of the column to check
        """
        # Check that column exists
        gfits = gammalib.GFits(filename)
        gtable = gfits[1]
        check = int(gtable.contains(colname))
        self.test_value(check, 1, 'Check column '+colname)

        # Return
        return


    def _check_normalization(self, filename, colnames):
        """
        Check that the probability values in columns identified by colnames are
        correctly normalized.

        Parameters
        ----------
        filename : str
            Event list file name
        colnames : list
            Names of the columns to check
        """
        # Check that the probability values are correctly normalized

        # Open fits table
        gfits  = gammalib.GFits(filename)
        gtable = gfits[1]
        ncols  = gtable.ncols()
        nevt   = gtable.nrows()
        srcs   = []

        # Loop over colnames
        for colname in colnames:

            # Get the correct column corresponding to colname from the table
            for icol in range(ncols):
                gcol = gtable[icol] 
                if gcol.name()==colname:
                    break

            # Read and save column content 
            src1 = []
            for ievt in range(nevt):
                src1.append(gcol[ievt])
            srcs.extend(src1)

        # Compute sum
        tot = sum(srcs)
        self.test_value(int(tot+0.5), nevt, 'Check that probability columns '
                                            'are normalized')

        # Return
        return

