#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csobsdef script.
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


# ============================== #
# Test class for csobsdef script #
# ============================== #
class Test(test):
    """
    Test class for csobsdef script

    This test class makes unit tests for the csobsdef script by using it
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
        self._pntdef_min = self._datadir + '/pntdef_minimal.dat'
        self._pntdef_max = self._datadir + '/pntdef_maximal.dat'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('csobsdef')

        # Append tests
        self.append(self._test_cmd, 'Test csobsdef on command line')
        self.append(self._test_python, 'Test csobsdef from Python')

        # Return
        return

    # Test csobsdef on command line
    def _test_cmd(self):
        """
        Test csobsdef on the command line
        """
        # Set script name
        csobsdef = self._script('csobsdef')

        # Setup csobsdef command
        cmd = csobsdef+' inpnt="'+self._pntdef_min+'"'+ \
                       ' outobs="csobsdef_cmd1.xml"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                       ' emin=0.1 emax=100.0 duration=1800.0 rad=5.0'+ \
                       ' logfile="csobsdef_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check observation definition file
        self._check_obsdef_file('csobsdef_cmd1.xml')

        # Setup csobsdef command
        cmd = csobsdef+' inpnt="obs_definition_that_does_not_exist.dat"'+ \
                       ' outobs="csobsdef_cmd1.xml"'+ \
                       ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                       ' emin=0.1 emax=100.0 duration=1800.0 rad=5.0'+ \
                       ' logfile="csobsdef_cmd2.log"'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test csobsdef from Python
    def _test_python(self):
        """
        Test csobsdef from Python
        """
        # Set-up csobsdef with minimal pointing definition file
        obsdef = cscripts.csobsdef()
        obsdef['inpnt']    = self._pntdef_min
        obsdef['outobs']   = 'csobsdef_py1.xml'
        obsdef['caldb']    = self._caldb
        obsdef['irf']      = self._irf
        obsdef['emin']     = 0.1
        obsdef['emax']     = 100.0
        obsdef['duration'] = 1800.0
        obsdef['rad']      = 5.0
        obsdef['logfile']  = 'csobsdef_py1.log'
        obsdef['chatter']  = 2

        # Run csobsdef script and save XML file
        obsdef.logFileOpen()   # Make sure we get a log file
        obsdef.run()
        obsdef.save()

        # Check model file
        self._check_obsdef_file('csobsdef_py1.xml')

        # Set-up csobsdef with maximal pointing definition file
        obsdef = cscripts.csobsdef()
        obsdef['inpnt']    = self._pntdef_max
        obsdef['outobs']   = 'csobsdef_py2.xml'
        obsdef['logfile']  = 'csobsdef_py2.log'
        obsdef['chatter']  = 3

        # Execute csobsdef script
        obsdef.execute()

        # Check model file
        self._check_obsdef_file('csobsdef_py2.xml')

        # Load CSV table
        csv = gammalib.GCsv(self._pntdef_min, ',')
 
        # Set-up csobsdef
        obsdef = cscripts.csobsdef()
        obsdef.pntdef(csv)
        obsdef['outobs']   = 'csobsdef_py3.xml'
        obsdef['caldb']    = self._caldb
        obsdef['irf']      = self._irf
        obsdef['emin']     = 0.1
        obsdef['emax']     = 100.0
        obsdef['duration'] = 1800.0
        obsdef['rad']      = 5.0
        obsdef['logfile']  = 'csobsdef_py3.log'
        obsdef['chatter']  = 4

        # Execute csobsdef script
        obsdef.execute()

        # Check model file
        self._check_obsdef_file('csobsdef_py3.xml')

        # Return
        return

    # Check observation definition file
    def _check_obsdef_file(self, filename):
        """
        Check observation definition file
        """
        # Load observation definition file
        obs = gammalib.GObservations(filename)

        # Check number of observations
        self.test_value(obs.size(), 3,
             'Check for 3 observations in observation definition file')

        # Loop over all observations
        for o in obs:
            self.test_value(o.instrument(), 'CTA',
                            'Check for "CTA" instrument')
            self.test_value(o.ontime(), 1800.0, 1.0e-6,
                            'Check for ontime')
            self.test_value(o.deadc(), 0.98, 1.0e-6,
                            'Check for deadtime correction')
            self.test_value(o.events().ebounds().emin().TeV(), 0.1, 1.0e-6,
                            'Check for minimum energy')
            self.test_value(o.events().ebounds().emax().TeV(), 100.0, 1.0e-6,
                            'Check for maximum energy')
            self.test_value(o.roi().radius(), 5.0, 1.0e-6,
                            'Check for ROI radius')

        # Return
        return
