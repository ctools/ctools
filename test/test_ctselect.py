#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctselect tool.
#
# Copyright (C) 2014-2016 Juergen Knoedlseder
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


# ============================ #
# Test class for ctselect tool #
# ============================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for ctselect tool.
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('ctselect')

        # Append tests
        self.append(self._test_cmd, 'Test ctselect on command line')
        self.append(self._test_python, 'Test ctselect from Python')

        # Return
        return

    # Test ctselect on command line
    def _test_cmd(self):
        """
        Test ctselect on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile('README.md'):
            ctselect = '../src/ctselect/ctselect'
        else:
            ctselect = 'ctselect'

        # Setup ctselect command
        cmd = ctselect+' inobs="data/crab_events.fits"'+ \
                       ' outobs="ctselect_cmd1.fits"'+ \
                       ' ra=83.63 dec=22.01 rad=3.0'+ \
                       ' tmin=0.0 tmax=1800.0'+ \
                       ' emin=0.1 emax=100.0'+ \
                       ' logfile="ctselect_cmd1.log" chatter=1'

        # Execute ctselect, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         'Successful ctselect execution on command line')

        # Check result file
        self._check_result_file('ctselect_cmd1.fits')

        # Setup ctselect command
        cmd = ctselect+' inobs="event_file_that_does_not_exist.fits"'+ \
                       ' outobs="ctselect_cmd2.fits"'+ \
                       ' ra=83.63 dec=22.01 rad=3.0'+ \
                       ' tmin=0.0 tmax=1800.0'+ \
                       ' emin=0.1 emax=100.0'+ \
                       ' logfile="ctselect_cmd2.log" chatter=1'

        # Execute ctselect, make sure we catch any exception
        try:
            rc = os.system(cmd+' >/dev/null 2>&1')
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         'Failure of ctselect execution on command line')

        # Return
        return

    # Test ctselect from Python
    def _test_python(self):
        """
        Test ctselect from Python.
        """
        # Set-up ctselect
        select = ctools.ctselect()
        select['inobs']   = 'data/crab_events.fits'
        select['outobs']  = 'ctselect_py1.fits'
        select['ra']      = 83.63
        select['dec']     = 22.01
        select['rad']     = 3
        select['tmin']    = 0
        select['tmax']    = 1800
        select['emin']    = 0.1
        select['emax']    = 100
        select['logfile'] = 'ctselect_py1.log'
        select['chatter'] = 2

        # Run ctselect tool
        select.logFileOpen()   # Make sure we get a log file
        select.run()
        select.save()

        # Check result file
        self._check_result_file('ctselect_py1.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file.
        """
        # Open result file
        result = gammalib.GCTAEventList(filename)

        # Check file
        self.test_value(result.size(), 6127, 'Check for 6127 events')
        self.test_value(result.roi().centre().dir().ra_deg(), 83.63, 1.0e-6,
                        'Check for ROI Right Ascension')
        self.test_value(result.roi().centre().dir().dec_deg(), 22.01, 1.0e-6,
                        'Check for ROI Declination')
        self.test_value(result.roi().radius(), 3.0, 1.0e-6,
                        'Check for ROI Radius')

        # Return
        return
