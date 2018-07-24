#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csobsinfo script.
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
import cscripts
from testing import test


# =============================== #
# Test class for csobsinfo script #
# =============================== #
class Test(test):
    """
    Test class for csobsinfo script

    This test class makes unit tests for the csobsinfo script by using it
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
        self._inobs = self._datadir + '/obs_unbinned.xml'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csobsinfo')

        # Append tests
        self.append(self._test_cmd, 'Test csobsinfo on command line')
        self.append(self._test_python, 'Test csobsinfo from Python')

        # Return
        return

    # Test csobsinfo on command line
    def _test_cmd(self):
        """
        Test csobsinfo on the command line
        """
        # Set script name
        csobsinfo = self._script('csobsinfo')

        # Setup csobsinfo command
        cmd = csobsinfo+' inobs="'+self._inobs+'"'+ \
                        ' outds9file="csobsinfo_cmd1.reg"'+ \
                        ' offset=yes ra=83.63 dec=22.01'+ \
                        ' logfile="csobsinfo_cmd1.log" chatter=1'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check DS9 file
        self._check_ds9_file('csobsinfo_cmd1.reg')

        # Setup csobsinfo command
        cmd = csobsinfo+' inobs="obs_definition_that_does_not_exist.xml"'+ \
                        ' outds9file="csobsinfo_cmd2.reg"'+ \
                        ' offset=yes ra=83.63 dec=22.01'+ \
                        ' logfile="csobsinfo_cmd2.log" debug=yes'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check csobsinfo --help
        self._check_help(csobsinfo)

        # Return
        return

    # Test csobsinfo from Python
    def _test_python(self):
        """
        Test csobsinfo from Python
        """
        # Set-up csobsinfo
        obsinfo = cscripts.csobsinfo()
        obsinfo['inobs']      = self._inobs
        obsinfo['outds9file'] = 'csobsinfo_py1.reg'
        obsinfo['offset']     = True
        obsinfo['ra']         = 83.63
        obsinfo['dec']        = 22.01
        obsinfo['logfile']    = 'csobsinfo_py1.log'
        obsinfo['chatter']    = 3

        # Run csobsinfo script and save XML file
        obsinfo.logFileOpen()   # Make sure we get a log file
        obsinfo.run()
        obsinfo.save()

        # Check DS9 file
        self._check_ds9_file('csobsinfo_py1.reg')

        # Check result access methods
        self.test_value(len(obsinfo.zeniths()), 1, 'Test number of zenith angles')
        self.test_value(obsinfo.zeniths()[0], 0.0, 'Test zenith angle')
        self.test_value(len(obsinfo.azimuths()), 1, 'Test number of azimuth angles')
        self.test_value(obsinfo.azimuths()[0], 0.0, 'Test azimuth angle')
        self.test_value(len(obsinfo.ras()), 1, 'Test number of Right Ascensions')
        self.test_value(obsinfo.ras()[0], 83.63, 'Test Right Ascension')
        self.test_value(len(obsinfo.decs()), 1, 'Test number of Declinations')
        self.test_value(obsinfo.decs()[0], 22.51, 'Test Declination')
        self.test_value(len(obsinfo.offsets()), 1, 'Test number of offsets')
        self.test_value(obsinfo.offsets()[0], 0.5, 'Test offset angle')
        self.test_value(obsinfo.ebounds().size(), 1, 'Test number of energy boundaries')
        self.test_value(obsinfo.ebounds().emin().TeV(), 1.0, 'Test minimum energy')
        self.test_value(obsinfo.ebounds().emax().TeV(), 100.0, 'Test maximum energy')

        # Return
        return

    # Check DS9 file
    def _check_ds9_file(self, filename):
        """
        Check DS9 file
        """
        # Open file   
        f = open(filename,'r')

        # Expect "fk5" in first line
        line = f.readline().strip('\n')
        self.test_value(line, 'fk5',
                         'Test for "fk5" in first line of DS9 file')

        # Expect "fk5" in first line
        line = f.readline().strip('\n')
        self.test_value(line[0:6], 'point(', 'Test for "point(" in second line of DS9 file')

        # Close file
        f.close()

        # Return
        return
