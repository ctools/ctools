#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csiactobs script.
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


# =============================== #
# Test class for csiactobs script #
# =============================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for csiactobs script.

    This test class makes unit tests for the csiactobs script by using it
    from the command line and from Python.
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
        self.name("csiactobs")

        # Append tests
        self.append(self._test_cmd, "Test csiactobs on command line")
        self.append(self._test_python, "Test csiactobs from Python")

        # Return
        return

    # Test csiactobs on command line
    def _test_cmd(self):
        """
        Test csiactobs on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile("README.md"):
            csiactobs = "../cscripts/csiactobs.py"
        else:
            csiactobs = "csiactobs"

        # Setup csiactobs command
        cmd = csiactobs+' datapath="iactdata"'+ \
                        ' prodname="unit-test"'+ \
                        ' infile="runlist_cmd1.dat"'+ \
                        ' bkgpars=1'+\
                        ' outobs="obs_iactobs_cmd1.xml"'+ \
                        ' outmodel="bgd_iactobs_cmd1.xml"'+ \
                        ' logfile="csiactobs_cmd1.log" chatter=1'

        # Execute csiactobs, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         "Successful csiactobs execution on command line")

        # Check observation definition XML file
        self._check_obsdef("obs_iactobs_cmd1.xml")

        # Check model definition XML file
        self._check_moddef("bgd_iactobs_cmd1.xml")

        # Setup csiactobs command
        cmd = csiactobs+' datapath="data_path_that_does_not_exist"'+ \
                        ' prodname="unit-test"'+ \
                        ' infile="runlist_cmd1.dat"'+ \
                        ' bkgpars=1'+\
                        ' outobs="obs_iactobs_cmd2.xml"'+ \
                        ' outmodel="bgd_iactobs_cmd2.xml"'+ \
                        ' logfile="csiactobs_cmd2.log" chatter=1'

        # Execute csiactobs, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         "Failure of csiactobs execution on command line")

        # Return
        return

    # Test csiactobs from Python
    def _test_python(self):
        """
        Test csiactobs from Python.
        """
        # Set-up csiactobs
        iactobs = cscripts.csiactobs()
        iactobs["datapath"] = "iactdata"
        iactobs["prodname"] = "unit-test"
        iactobs["infile"]   = "runlist_py1.dat"
        iactobs["bkgpars"]  = 1
        iactobs["outobs"]   = "obs_iactobs_py1.xml"
        iactobs["outmodel"] = "bgd_iactobs_py1.xml"
        iactobs["logfile"]  = "csiactobs_py1.log"
        iactobs["chatter"]  = 2

        # Run csiactobs script and save run list
        iactobs.logFileOpen()   # Make sure we get a log file
        iactobs.run()
        iactobs.save()

        # Check observation definition XML file
        self._check_obsdef("obs_iactobs_py1.xml")
        
        # Check model definition XML file
        self._check_moddef("bgd_iactobs_py1.xml")

        # Return
        return

    # Check observation definition XML file
    def _check_obsdef(self, filename):
        """
        Check observation definition XML file.
        """
        # Load observation definition XML file
        obs = gammalib.GObservations(filename)

        # Get response
        rsp = obs[0].response()

        # Check number of observations
        self.test_value(obs.size(), 1,
                        'Check for single observation in XML file')
        self.test_assert(obs[0].eventfile().file() == 'crab_events.fits.gz',
                         'Check event file name')
        self.test_assert(obs[0].eventfile().extname() == 'EVENTS',
                         'Check event extension name')
        self.test_assert(rsp.aeff().filename().file() ==
                         'irf_file.fits.gz',
                         'Check effective area name')
        self.test_assert(rsp.aeff().filename().extname() ==
                         'EFFECTIVE AREA',
                         'Check effective area extension name')
        self.test_assert(rsp.psf().filename().file() ==
                         'irf_file.fits.gz',
                         'Check point spread function name')
        self.test_assert(rsp.psf().filename().extname() ==
                         'POINT SPREAD FUNCTION',
                         'Check point spread function extension name')
        self.test_assert(rsp.edisp().filename().file() ==
                         'irf_file.fits.gz',
                         'Check energy dispersion name')
        self.test_assert(rsp.edisp().filename().extname() ==
                         'ENERGY DISPERSION',
                         'Check energy dispersion extension name')
        self.test_assert(rsp.background().filename().file() ==
                         'irf_file.fits.gz',
                         'Check background name')
        self.test_assert(rsp.background().filename().extname() ==
                         'BACKGROUND',
                         'Check background extension name')
        
        # Return
        return

    # Check model XML file
    def _check_moddef(self, filename):
        """
        Check model definition XML file.
        """
        # Load model definition XML file
        models = gammalib.GModels(filename)

        # Check number of models
        self.test_value(models.size(), 1,
                        'Check for single model in XML file')

        # Return
        return
