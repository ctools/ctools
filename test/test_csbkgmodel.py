#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csbkgmodel script.
#
# Copyright (C) 2019 Juergen Knoedlseder
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
# Test class for csbkgmodel script #
# ================================ #
class Test(test):
    """
    Test class for csbkgmodel script

    This test class makes unit tests for the csbkgmodel script by using it
    from the command line and from Python.
    """

    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        test.__init__(self)

        # Set test datasets and parameters
        self._myevents1 = self._datadir + '/crab_offaxis1.fits'
        self._myevents2 = self._datadir + '/crab_offaxis2.fits'
        self._lookup    = self._datadir + '/hess_bkg_lookup.fits'

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('csbkgmodel')

        # Append tests
        self.append(self._test_cmd, 'Test csbkgmodel on command line')
        self.append(self._test_python, 'Test csbkgmodel from Python')
        self.append(self._test_pickeling, 'Test csbkgmodel pickeling')

        # Return
        return

    # Test csbkgmodel on command line
    def _test_cmd(self):
        """
        Test csbkgmodel on the command line.
        """
        # Set script name
        csbkgmodel = self._script('csbkgmodel')

        # Setup csbkgmodel command
        cmd = csbkgmodel+' inobs="'+self._events+'"'+ \
                         ' caldb="'+self._caldb+'" irf="'+self._irf+'"' \
                         ' instrument=CTA spatial=GAUSS gradient=yes'+ \
                         ' spectral=NODES ebinalg=LOG emin=1.0 emax=100.0'+ \
                         ' enumbins=8 runwise=yes rad=2.0'+ \
                         ' outmodel="csbkgmodel_cmd1.xml"'+ \
                         ' logfile="csbkgmodel_cmd1.log" chatter=4'

        # Check if execution was successful
        self.test_assert(self._execute(cmd) == 0,
             'Check successful execution from command line')

        # Check background model
        self._check_bkg_model('csbkgmodel_cmd1.xml')

        # Setup csbkgmodel command
        cmd = csbkgmodel+' inobs="events_that_do_not_exist.fits"'+ \
                         ' caldb="'+self._caldb+'" irf="'+self._irf+'"' \
                         ' instrument=CTA spatial=GAUSS gradient=yes'+ \
                         ' spectral=NODES ebinalg=LOG emin=1.0 emax=100.0'+ \
                         ' enumbins=8 runwise=yes rad=2.0'+ \
                         ' outmodel="csbkgmodel_cmd2.xml"' + \
                         ' logfile="csbkgmodel_cmd2.log" chatter=1'

        # Check if execution failed
        self.test_assert(self._execute(cmd, success=False) != 0,
             'Check invalid input file when executed from command line')

        # Check csbkgmodel --help
        self._check_help(csbkgmodel)

        # Return
        return

    # Test csbkgmodel from Python
    def _test_python(self):
        """
        Test csbkgmodel from Python
        """
        # Set-up csbkgmodel
        bkgmodel = cscripts.csbkgmodel()
        bkgmodel['inobs']      = self._events
        bkgmodel['caldb']      = self._caldb
        bkgmodel['irf']        = self._irf
        bkgmodel['instrument'] = 'CTA'
        bkgmodel['spatial']    = 'GAUSS'
        bkgmodel['gradient']   = True
        bkgmodel['spectral']   = 'NODES'
        bkgmodel['ebinalg']    = 'LOG'
        bkgmodel['emin']       = 1.0
        bkgmodel['emax']       = 100.0
        bkgmodel['enumbins']   = 8
        bkgmodel['runwise']    = True
        bkgmodel['rad']        = 2.0
        bkgmodel['chatter']    = 2
        bkgmodel['outmodel']   = 'csbkgmodel_py1.xml'
        bkgmodel['logfile']    = 'csbkgmodel_py1.log'

        # Run csbkgmodel script and save background model
        bkgmodel.logFileOpen()   # Make sure we get a log file
        bkgmodel.run()
        bkgmodel.save()

        # Check background model
        self._check_bkg_model('csbkgmodel_py1.xml')

        # Now test without gradient, power law and not runwise
        bkgmodel = cscripts.csbkgmodel()
        bkgmodel['inobs']      = self._events
        bkgmodel['caldb']      = self._caldb
        bkgmodel['irf']        = self._irf
        bkgmodel['instrument'] = 'CTA'
        bkgmodel['spatial']    = 'GAUSS'
        bkgmodel['gradient']   = False
        bkgmodel['spectral']   = 'PLAW'
        bkgmodel['ebinalg']    = 'LOG'
        bkgmodel['emin']       = 1.0
        bkgmodel['emax']       = 100.0
        bkgmodel['enumbins']   = 8
        bkgmodel['runwise']    = False
        bkgmodel['rad']        = 2.0
        bkgmodel['chatter']    = 3
        bkgmodel['outmodel']   = 'csbkgmodel_py2.xml'
        bkgmodel['logfile']    = 'csbkgmodel_py2.log'

        # Execute csbkgmodel script
        bkgmodel.logFileOpen()   # Make sure we get a log file
        bkgmodel.execute()

        # Check background model
        self._check_bkg_model('csbkgmodel_py2.xml')

        # Now test AEFF model
        bkgmodel = cscripts.csbkgmodel()
        bkgmodel['inobs']      = self._events
        bkgmodel['caldb']      = self._caldb
        bkgmodel['irf']        = self._irf
        bkgmodel['instrument'] = 'CTA'
        bkgmodel['spatial']    = 'AEFF'
        bkgmodel['gradient']   = False
        bkgmodel['spectral']   = 'PLAW'
        bkgmodel['ebinalg']    = 'LOG'
        bkgmodel['emin']       = 1.0
        bkgmodel['emax']       = 100.0
        bkgmodel['enumbins']   = 8
        bkgmodel['runwise']    = False
        bkgmodel['rad']        = 2.0
        bkgmodel['chatter']    = 4
        bkgmodel['outmodel']   = 'csbkgmodel_py3.xml'
        bkgmodel['logfile']    = 'csbkgmodel_py3.log'

        # Execute csbkgmodel script
        bkgmodel.logFileOpen()   # Make sure we get a log file
        bkgmodel.execute()

        # Check background model
        self._check_bkg_model('csbkgmodel_py3.xml')

        # Now test IRF model
        bkgmodel = cscripts.csbkgmodel()
        bkgmodel['inobs']      = self._events
        bkgmodel['caldb']      = self._caldb
        bkgmodel['irf']        = self._irf
        bkgmodel['instrument'] = 'CTA'
        bkgmodel['spatial']    = 'IRF'
        bkgmodel['gradient']   = False
        bkgmodel['spectral']   = 'PLAW'
        bkgmodel['ebinalg']    = 'LOG'
        bkgmodel['emin']       = 1.0
        bkgmodel['emax']       = 100.0
        bkgmodel['enumbins']   = 8
        bkgmodel['runwise']    = False
        bkgmodel['rad']        = 2.0
        bkgmodel['chatter']    = 4
        bkgmodel['outmodel']   = 'csbkgmodel_py4.xml'
        bkgmodel['logfile']    = 'csbkgmodel_py4.log'

        # Execute csbkgmodel script
        bkgmodel.logFileOpen()   # Make sure we get a log file
        bkgmodel.execute()

        # Check background model
        self._check_bkg_model('csbkgmodel_py4.xml')

        # Test with multiple input observations
        obs = gammalib.GObservations()
        for s, events in enumerate([self._myevents1, self._myevents2]):
            run = gammalib.GCTAObservation(events)
            run.id(str(s + 1))
            run.response(self._irf, gammalib.GCaldb('cta', self._caldb))
            obs.append(run)

        # Set-up csbkgmodel
        bkgmodel = cscripts.csbkgmodel(obs)
        bkgmodel['instrument'] = 'CTA'
        bkgmodel['spatial']    = 'GAUSS'
        bkgmodel['gradient']   = True
        bkgmodel['spectral']   = 'NODES'
        bkgmodel['ebinalg']    = 'POW'
        bkgmodel['emin']       = 1.0
        bkgmodel['emax']       = 100.0
        bkgmodel['enumbins']   = 8
        bkgmodel['ebingamma']  = 1.1
        bkgmodel['runwise']    = True
        bkgmodel['rad']        = 2.0
        bkgmodel['chatter']    = 2
        bkgmodel['outmodel']   = 'csbkgmodel_py5.xml'
        bkgmodel['logfile']    = 'csbkgmodel_py5.log'

        # Execute csbkgmodel script
        bkgmodel.logFileOpen()   # Make sure we get a log file
        bkgmodel.execute()

        # Check background model
        self._check_bkg_model('csbkgmodel_py5.xml', nmodels=2)

        # Test GAUSS(E) spatial model
        bkgmodel = cscripts.csbkgmodel()
        bkgmodel['inobs']      = self._events
        bkgmodel['caldb']      = self._caldb
        bkgmodel['irf']        = self._irf
        bkgmodel['instrument'] = 'CTA'
        bkgmodel['spatial']    = 'GAUSS(E)'
        bkgmodel['snumbins']   = 2
        bkgmodel['smin']       = 1.0
        bkgmodel['smax']       = 10.0
        bkgmodel['gradient']   = True
        bkgmodel['spectral']   = 'NODES'
        bkgmodel['ebinalg']    = 'LOG'
        bkgmodel['emin']       = 1.0
        bkgmodel['emax']       = 100.0
        bkgmodel['enumbins']   = 8
        bkgmodel['runwise']    = True
        bkgmodel['rad']        = 2.0
        bkgmodel['chatter']    = 2
        bkgmodel['outmodel']   = 'csbkgmodel_py6.xml'
        bkgmodel['logfile']    = 'csbkgmodel_py6.log'

        # Run csbkgmodel script and save background model
        bkgmodel.logFileOpen()   # Make sure we get a log file
        bkgmodel.execute()

        # Check background model
        self._check_bkg_model('csbkgmodel_py6.xml')

        # Test LOOKUP spatial model
        bkgmodel = cscripts.csbkgmodel()
        bkgmodel['inobs']      = self._events
        bkgmodel['caldb']      = self._caldb
        bkgmodel['irf']        = self._irf
        bkgmodel['instrument'] = 'CTA'
        bkgmodel['spatial']    = 'LOOKUP'
        bkgmodel['slufile']    = self._lookup
        bkgmodel['gradient']   = True
        bkgmodel['spectral']   = 'NODES'
        bkgmodel['ebinalg']    = 'LOG'
        bkgmodel['emin']       = 1.0
        bkgmodel['emax']       = 100.0
        bkgmodel['enumbins']   = 8
        bkgmodel['runwise']    = True
        bkgmodel['rad']        = 2.0
        bkgmodel['chatter']    = 2
        bkgmodel['outmodel']   = 'csbkgmodel_py7.xml'
        bkgmodel['logfile']    = 'csbkgmodel_py7.log'

        # Run csbkgmodel script and save background model
        bkgmodel.logFileOpen()   # Make sure we get a log file
        bkgmodel.execute()

        # Check background model
        self._check_bkg_model('csbkgmodel_py7.xml')

        # Test PROFILE spatial model
        bkgmodel = cscripts.csbkgmodel()
        bkgmodel['inobs']      = self._events
        bkgmodel['caldb']      = self._caldb
        bkgmodel['irf']        = self._irf
        bkgmodel['instrument'] = 'CTA'
        bkgmodel['spatial']    = 'PROFILE'
        bkgmodel['gradient']   = True
        bkgmodel['spectral']   = 'NODES'
        bkgmodel['ebinalg']    = 'LOG'
        bkgmodel['emin']       = 1.0
        bkgmodel['emax']       = 100.0
        bkgmodel['enumbins']   = 8
        bkgmodel['runwise']    = True
        bkgmodel['rad']        = 2.0
        bkgmodel['chatter']    = 2
        bkgmodel['outmodel']   = 'csbkgmodel_py8.xml'
        bkgmodel['logfile']    = 'csbkgmodel_py8.log'

        # Run csbkgmodel script and save background model
        bkgmodel.logFileOpen()   # Make sure we get a log file
        bkgmodel.execute()

        # Check background model
        self._check_bkg_model('csbkgmodel_py8.xml')

        # Test POLYNOM spatial model
        bkgmodel = cscripts.csbkgmodel()
        bkgmodel['inobs']      = self._events
        bkgmodel['caldb']      = self._caldb
        bkgmodel['irf']        = self._irf
        bkgmodel['instrument'] = 'CTA'
        bkgmodel['spatial']    = 'POLYNOM'
        bkgmodel['gradient']   = True
        bkgmodel['spectral']   = 'NODES'
        bkgmodel['ebinalg']    = 'LOG'
        bkgmodel['emin']       = 1.0
        bkgmodel['emax']       = 100.0
        bkgmodel['enumbins']   = 8
        bkgmodel['runwise']    = True
        bkgmodel['rad']        = 2.0
        bkgmodel['chatter']    = 2
        bkgmodel['outmodel']   = 'csbkgmodel_py9.xml'
        bkgmodel['logfile']    = 'csbkgmodel_py9.log'

        # Run csbkgmodel script and save background model
        bkgmodel.logFileOpen()   # Make sure we get a log file
        bkgmodel.execute()

        # Check background model
        self._check_bkg_model('csbkgmodel_py9.xml')

        # Return
        return

    # Test csbkgmodel pickeling
    def _test_pickeling(self):
        """
        Test csbkgmodel pickeling
        """
        # Perform pickeling tests of empty class
        self._pickeling(cscripts.csbkgmodel())

        # Set-up csbkgmodel
        bkgmodel = cscripts.csbkgmodel()
        bkgmodel['inobs']      = self._events
        bkgmodel['caldb']      = self._caldb
        bkgmodel['irf']        = self._irf
        bkgmodel['instrument'] = 'CTA'
        bkgmodel['spatial']    = 'GAUSS'
        bkgmodel['gradient']   = True
        bkgmodel['spectral']   = 'NODES'
        bkgmodel['ebinalg']    = 'LOG'
        bkgmodel['emin']       = 1.0
        bkgmodel['emax']       = 100.0
        bkgmodel['enumbins']   = 8
        bkgmodel['runwise']    = True
        bkgmodel['rad']        = 2.0
        bkgmodel['chatter']    = 2
        bkgmodel['outmodel']   = 'csbkgmodel_py1_pickle.xml'
        bkgmodel['logfile']    = 'csbkgmodel_py1_pickle.log'

        # Perform pickeling tests of filled class
        obj = self._pickeling(bkgmodel)

        # Run csbkgmodel script and save background model
        obj.logFileOpen()   # Make sure we get a log file
        obj.run()
        obj.save()

        # Check background model
        self._check_bkg_model('csbkgmodel_py1_pickle.xml')

        # Return
        return

    # Check background model
    def _check_bkg_model(self, filename, nmodels=1):
        """
        Check background model
        """
        # Load model
        models = gammalib.GModels(filename)

        # Check model size
        self.test_value(models.size(), nmodels,
                        'Check for %d model(s) in XML file' % nmodels)

        # Return
        return
