#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csobsselect script
#
# Copyright (C) 2016 Juergen Knoedlseder
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


# ================================= #
# Test class for csobsselect script #
# ================================= #
class Test(test):
    """
    Test class for csobsselect script

    This test class makes unit tests for the csobsselect script by using it
    from the command line and from Python.
    """
    
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        test.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('csobsselect')

        # Append tests
        self.append(self._test_cmd, 'Test csobsselect on command line')
        self.append(self._test_python, 'Test csobsselect from Python')

        # Return
        return

    # Test csobsselect on command line
    def _test_cmd(self):
        """
        Test csobsselect on the command line
        """
        # Set script name
        csobsselect = self._script('csobsselect')

        # Setup csobsselect command
        cmd = csobsselect+' inobs="'+self._events+'"'+ \
                          ' outobs="csobsselect_cmd1.xml"'+ \
                          ' pntselect="CIRCLE" coordsys="CEL" ra=83.63'+ \
                          ' dec=22.01 rad=5.0'+ \
                          ' logfile="csobsselect_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check observation definition file
        self._check_observation_file('csobsselect_cmd1.xml', 1)

        # Setup csobsselect command
        cmd = csobsselect+' inobs="observation_that_does_not_exist.xml"'+ \
                          ' outobs="csobsselect_cmd2.xml"'+ \
                          ' pntselect="CIRCLE" coordsys="CEL" ra=83.63'+ \
                          ' dec=22.01 rad=5.0'+ \
                          ' logfile="csobsselect_cmd2.log" chatter=2'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csobsselect from Python
    def _test_python(self):
        """
        Test csobsselect from Python
        """
        # Allocate empty csobsselect script
        obsselect = cscripts.csobsselect()

        # Check that saving does not nothing
        obsselect['outobs']  = 'csobsselect_py0.xml'
        obsselect['logfile'] = 'csobsselect_py0.log'
        obsselect.logFileOpen()
        obsselect.save()

        # Check than an empty observation definition file has been created
        self._check_observation_file('csobsselect_py0.xml', 0)

        # Check that clearing does not lead to an exception or segfault
        obsselect.clear()

        # Now set csobsselect parameters
        obsselect['inobs']     = self._events
        obsselect['outobs']    = 'csobsselect_py1.xml'
        obsselect['pntselect'] = 'CIRCLE'
        obsselect['coordsys']  = 'CEL'
        obsselect['ra']        = 83.63
        obsselect['dec']       = 22.01
        obsselect['rad']       = 5.0
        obsselect['logfile']   = 'csobsselect_py1.log'
        obsselect['chatter']   = 2

        # Run csobsselect script and save observations
        obsselect.logFileOpen()   # Make sure we get a log file
        obsselect.run()
        obsselect.save()

        # Check model file
        self._check_observation_file('csobsselect_py1.xml', 1)

        # Execute csobsselect script again, now in Galactic coordinates
        obsselect['outobs']    = 'csobsselect_py2.xml'
        obsselect['pntselect'] = 'CIRCLE'
        obsselect['coordsys']  = 'GAL'
        obsselect['glon']      = 184.56
        obsselect['glat']      = -5.79
        obsselect['rad']       = 5.0
        obsselect['logfile']   = 'csobsselect_py2.log'
        obsselect['chatter']   = 3
        obsselect.logFileOpen()  # Needed to get a new log file
        obsselect.execute()

        # Check model file
        self._check_observation_file('csobsselect_py2.xml', 1)

        # Execute csobsselect script again, now a box in Galactic coordinates
        obsselect['outobs']    = 'csobsselect_py3.xml'
        obsselect['pntselect'] = 'BOX'
        obsselect['coordsys']  = 'GAL'
        obsselect['glon']      = 184.56
        obsselect['glat']      = -5.79
        obsselect['width']     = 5.0
        obsselect['height']    = 5.0
        obsselect['logfile']   = 'csobsselect_py3.log'
        obsselect['chatter']   = 4
        obsselect.logFileOpen()  # Needed to get a new log file
        obsselect.execute()

        # Check model file
        self._check_observation_file('csobsselect_py3.xml', 1)

        # Execute csobsselect script again, now a box in celestial coordinates
        obsselect['outobs']    = 'csobsselect_py4.xml'
        obsselect['pntselect'] = 'BOX'
        obsselect['coordsys']  = 'CEL'
        obsselect['ra']        = 83.63
        obsselect['dec']       = 22.01
        obsselect['width']     = 5.0
        obsselect['height']    = 5.0
        obsselect['logfile']   = 'csobsselect_py4.log'
        obsselect['chatter']   = 4
        obsselect.logFileOpen()  # Needed to get a new log file
        obsselect.execute()

        # Check model file
        self._check_observation_file('csobsselect_py4.xml', 1)

        # Return
        return

    # Check observation definition file
    def _check_observation_file(self, filename, number):
        """
        Check observation definition file
        """
        # Open observations
        obs = gammalib.GObservations(filename)

        # Check number of models
        self.test_value(obs.size(), number, 'Check for '+str(number)+
                        ' observations in XML file')
        
        # Return
        return
