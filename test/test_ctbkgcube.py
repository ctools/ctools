#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the ctbkgcube tool.
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
import gammalib
import ctools
from testing import test


# ============================= #
# Test class for ctbkgcube tool #
# ============================= #
class Test(test):
    """
    Test class for ctbkgcube tool.
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('ctbkgcube')

        # Append tests
        self.append(self._test_cmd, 'Test ctbkgcube on command line')
        self.append(self._test_python, 'Test ctbkgcube from Python')

        # Return
        return

    # Test ctbkgcube on command line
    def _test_cmd(self):
        """
        Test ctbkgcube on the command line.
        """
        # Set tool name
        ctbkgcube = self._tool('ctbkgcube')

        # Setup ctbkgcube command
        cmd = ctbkgcube+' inobs="'+self._events+'"'+ \
                        ' inmodel="'+self._model+'"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctbkgcube_cmd1.fits"'+ \
                        ' outmodel="ctbkgcube_cmd1.xml"'+ \
                        ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=10 nypix=10 binsz=0.4'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' logfile="ctbkgcube_cmd1.log" chatter=1'

        # Check if execution of wrong command fails
        self.test_assert(self._execute('command_that_does_not_exist') != 0,
             'Self test of test script')

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check result file
        self._check_result_file('ctbkgcube_cmd1.fits')

        # Setup ctbkgcube command
        cmd = ctbkgcube+' inobs="event_file_that_does_not_exist.fits"'+ \
                        ' inmodel="'+self._model+'"'+ \
                        ' incube="NONE"'+ \
                        ' outcube="ctbkgcube_cmd2.fits"'+ \
                        ' outmodel="ctbkgcube_cmd2.xml"'+ \
                        ' caldb="'+self._caldb+'" irf="'+self._irf+'"'+ \
                        ' ebinalg="LOG" emin=0.1 emax=100.0 enumbins=20'+ \
                        ' nxpix=10 nypix=10 binsz=0.4'+ \
                        ' coordsys="CEL" proj="CAR" xref=83.63 yref=22.01'+ \
                        ' logfile="ctbkgcube_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd) != 0,
             'Check invalid input file when executed from command line')

        # Return
        return

    # Test ctbkgcube from Python
    def _test_python(self):
        """
        Test ctbkgcube from Python
        """
        # Set-up ctbkgcube
        bkgcube = ctools.ctbkgcube()
        bkgcube['inobs']    = self._events
        bkgcube['inmodel']  = self._model
        bkgcube['incube']   = 'NONE'
        bkgcube['caldb']    = self._caldb
        bkgcube['irf']      = self._irf
        bkgcube['ebinalg']  = 'LOG'
        bkgcube['emin']     = 0.1
        bkgcube['emax']     = 100
        bkgcube['enumbins'] = 20
        bkgcube['nxpix']    = 10
        bkgcube['nypix']    = 10
        bkgcube['binsz']    = 0.4
        bkgcube['coordsys'] = 'CEL'
        bkgcube['proj']     = 'CAR'
        bkgcube['xref']     = 83.63
        bkgcube['yref']     = 22.01
        bkgcube['outcube']  = 'ctbkgcube_py1.fits'
        bkgcube['outmodel'] = 'ctbkgcube_py1.xml'
        bkgcube['logfile']  = 'ctbkgcube_py1.log'
        bkgcube['chatter']  = 2

        # Run ctbkgcube tool
        bkgcube.logFileOpen()   # Make sure we get a log file
        bkgcube.run()
        bkgcube.save()

        # Check result file
        self._check_result_file('ctbkgcube_py1.fits')

        # Set-up ctbkgcube from counts cube
        bkgcube = ctools.ctbkgcube()
        bkgcube['inobs']     = self._events
        bkgcube['inmodel']   = self._model
        bkgcube['incube']    = self._cntcube
        bkgcube['caldb']     = self._caldb
        bkgcube['irf']       = self._irf
        bkgcube['outcube']   = 'ctbkgcube_py2.fits'
        bkgcube['outmodel']  = 'ctbkgcube_py2.xml'
        bkgcube['logfile']   = 'ctbkgcube_py2.log'
        bkgcube['chatter']   = 3
        bkgcube['publish']   = True
        bkgcube['addbounds'] = True

        # Execute ctbkgcube tool
        bkgcube.logFileOpen()   # Make sure we get a log file
        bkgcube.execute()

        # Check result file
        self._check_result_file('ctbkgcube_py2.fits')

        # Copy ctbkgcube tool and execute copied tool
        cpy_bkgcube = bkgcube
        cpy_bkgcube['outcube']  = 'ctbkgcube_py3.fits'
        cpy_bkgcube['outmodel'] = 'ctbkgcube_py3.xml'
        cpy_bkgcube['logfile']  = 'ctbkgcube_py3.log'
        cpy_bkgcube['chatter']  = 4
        cpy_bkgcube.execute()

        # Check result file
        self._check_result_file('ctbkgcube_py3.fits')

        # Clear ctbkgcube tool
        bkgcube.clear()

        # TODO: Do some test after clearing of ctbkgcube tool

        # Prepare observation container
        obs1 = gammalib.GCTAObservation(self._events)
        obs2 = gammalib.GCTAObservation(self._cntcube)
        obs3 = gammalib.GLATObservation()
        obs1.id('0001')
        obs2.id('0002')
        obs3.id('0001')
        obs3.events(gammalib.GLATEventList())
        obs = gammalib.GObservations()
        obs.append(obs1)
        obs.append(obs2)
        obs.append(obs3)
        obs.models(gammalib.GModels(self._model))

        # Set-up ctbkgcube from observation container
        bkgcube = ctools.ctbkgcube(obs)
        #bkgcube['inobs']    = self._events
        #bkgcube['inmodel']  = self._model
        bkgcube['incube']   = 'NONE'
        bkgcube['caldb']    = self._caldb
        bkgcube['irf']      = self._irf
        bkgcube['ebinalg']  = 'LOG'
        bkgcube['emin']     = 0.1
        bkgcube['emax']     = 100
        bkgcube['enumbins'] = 20
        bkgcube['nxpix']    = 10
        bkgcube['nypix']    = 10
        bkgcube['binsz']    = 0.4
        bkgcube['coordsys'] = 'CEL'
        bkgcube['proj']     = 'CAR'
        bkgcube['xref']     = 83.63
        bkgcube['yref']     = 22.01
        bkgcube['outcube']  = 'ctbkgcube_py4.fits'
        bkgcube['outmodel'] = 'ctbkgcube_py4.xml'
        bkgcube['logfile']  = 'ctbkgcube_py4.log'
        bkgcube['chatter']  = 4

        # Execute ctbkgcube tool
        bkgcube.logFileOpen()   # Make sure we get a log file
        bkgcube.execute()

        # Check result file
        self._check_result_file('ctbkgcube_py4.fits')

        # Return
        return

    # Check result file
    def _check_result_file(self, filename):
        """
        Check result file
        """
        # Open result file
        result = gammalib.GCTACubeBackground(filename)

        # Check dimensions
        self.test_value(len(result.energies()), 21, 'Check for 21 energy values')

        # Return
        return
