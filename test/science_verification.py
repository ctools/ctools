#! /usr/bin/env python
# ==========================================================================
# This script performs the ctools science verification. It creates and
# analyses the pull distributions for a variety of spectral and spatial
# models. Test are generally done in unbinned mode, but also a stacked
# analysis test is included. At the end the script produces a JUnit
# compliant science verification report.
#
# Usage:
#   ./science_verification.py
#
# --------------------------------------------------------------------------
#
# Copyright (C) 2015-2017 Juergen Knoedlseder
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
import csv
import math
import sys
import gammalib
import ctools
import cscripts


# ========================== #
# Generate pull distribution #
# ========================== #
def generate_pull_distribution(model, obs='NONE', trials=100, \
                               caldb='prod2', irf='South_50h', \
                               deadc=0.98, edisp=False, \
                               ra=83.63, dec=22.01, rad=5.0, \
                               emin=0.1, emax=100.0, enumbins=0, \
                               duration=1800.0, \
                               npix=200, binsz=0.02, \
                               coordsys='CEL', proj='TAN', \
                               pattern='single', offset=1.5, \
                               debug=False, chatter=2):
    """
    Generates pull distribution for a given model

    Parameters
    ----------
    model : str
        Model XML filename (without .xml extension)
    obs : str, optional
        Input observation definition XML filename
    trials : int, optional
        Number of trials
    caldb : str, optional
        Calibration database
    irf : str, optional
        Name of instrument response function
    deadc : float, optional
        Deadtime correction factor
    edisp : bool, optional
        Use energy dispersion?
    ra : float, optional
        Right Ascension of pointing (deg)
    dec : float, optional
        Declination of pointing (deg)
    rad : float, optional
        Simulation radius (deg)
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    enumbins : int, optional
        Number of energy bins (0 for unbinned analysis)
    duration : float, optional
        Observation duration (sec)
    npix : int, optional
        Number of pixels
    binsz : float, optional
        Pixel size (deg/pixel)
    coordsys : str, optional
        Coordinate system (CEL or GAL)
    proj : str, optional
        Sky projection
    pattern : str, optional
        Pointing pattern
    offset : float, optional
        Offset between pointing directions
    debug : bool, optional
        Enable debugging?
    chatter : int, optional
        Chatter level

    Returns
    -------
    outfile : str
        Name of pull distribution output file
    """
    # Derive parameters
    head, tail = os.path.split(model)
    inmodel    = model + '.xml'
    outfile    = 'cspull_' + tail + '.dat'

    # Setup pull distribution generation
    pull = cscripts.cspull()
    pull['inobs']    = obs
    pull['inmodel']  = inmodel
    pull['outfile']  = outfile
    pull['caldb']    = caldb
    pull['irf']      = irf
    pull['edisp']    = edisp
    pull['ra']       = ra
    pull['dec']      = dec
    pull['rad']      = rad
    pull['emin']     = emin
    pull['emax']     = emax
    pull['tmin']     = 0.0
    pull['tmax']     = duration
    pull['enumbins'] = enumbins
    pull['npix']     = npix
    pull['binsz']    = binsz
    pull['coordsys'] = coordsys
    pull['proj']     = proj
    pull['deadc']    = deadc
    pull['rad']      = rad
    pull['pattern']  = pattern
    pull['offset']   = offset
    pull['ntrials']  = trials
    pull['debug']    = debug
    pull['chatter']  = chatter

    # Generate pull distributions
    pull.run()

    # Return
    return outfile


# ========================= #
# Analyse pull distribution #
# ========================= #
def analyse_pull_distribution(filename):
    """
    Compute mean and standard deviation of pull distribution

    Parameters
    ----------
    filename : str
        Pull distribution ASCII file to analyse

    Returns
    -------
    results : dict
        Result dictionary
    """
    # Initialise column names, means and standard deviations
    colnames = []
    means    = []
    stds     = []

    # Open reader
    reader = csv.reader(open(filename, 'r'), delimiter=',')

    # Read rows
    first   = True
    index   = -1
    samples = 0.0
    for row in reader:

        # Get column names if first row
        if first:
            for element in row:
                colnames.append(element)
                means.append(0.0)
                stds.append(0.0)

        # Handle data rows
        else:
            for i, element in enumerate(row):
                means[i] += float(element)
                stds[i]  += float(element)*float(element)
            samples += 1.0

        # Flag that first row has been passed
        first = False

    # Compute mean and standard deviations
    for i in range(len(stds)):
        std       = math.sqrt(stds[i]/samples -
                              means[i]*means[i]/(samples*samples))
        stds[i]   = std
        means[i] /= samples

    # Setup results
    results = {}
    for i in range(len(colnames)):
        results[colnames[i]] = {'mean': means[i], 'std':  stds[i]}

    # Return results
    return results


# =================================== #
# Test class for science verification #
# =================================== #
class sciver(gammalib.GPythonTestSuite):
    """
    Test class for science verification
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Initialise results
        self.results = None

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('Science Verification')

        # Append background model test
        self.append(self.bgd, 'Test background model')

        # Append spectral tests
        self.append(self.spec_plaw, 'Test power law model')
        self.append(self.spec_plaw_edisp, 'Test power law model with energy dispersion')
        self.append(self.spec_plaw_stacked, 'Test power law model with stacked analysis')
        self.append(self.spec_plaw2, 'Test power law 2 model')
        self.append(self.spec_smoothbplaw, 'Test smoothly broken power law model')
        self.append(self.spec_eplaw, 'Test exponentially cut off power law model')
        self.append(self.spec_supeplaw, 'Test super exponentially cut off power law model')
        self.append(self.spec_logparabola, 'Test log parabola model')
        self.append(self.spec_gauss, 'Test Gaussian model')
        self.append(self.spec_filefct, 'Test file function model')
        self.append(self.spec_nodes, 'Test nodes model')

        # Append spatial tests
        self.append(self.spat_ptsrc, 'Test point source model')
        self.append(self.spat_rdisk, 'Test radial disk model')
        self.append(self.spat_rgauss, 'Test radial Gaussian model')
        self.append(self.spat_rshell, 'Test radial shell model')
        self.append(self.spat_edisk, 'Test elliptical disk model')
        self.append(self.spat_egauss, 'Test elliptical Gaussian model')
        self.append(self.spat_const, 'Test diffuse isotropic model')
        self.append(self.spat_map, 'Test diffuse map model')
        self.append(self.spat_map_roi, 'Test diffuse map model (small ROI)')
        self.append(self.spat_map_nn, 'Test diffuse map model (not normalized and scaled)')
        self.append(self.spat_cube, 'Test diffuse cube model')

        # Return
        return

    # Generate and analyse pull distributions
    def pull(self, model, obs='NONE', trials=100, duration=1800.0,
             ra=83.63, dec=22.01, rad=5.0,
             emin=0.1, emax=100.0, enumbins=0,
             edisp=False, debug=False):
        """
        Generate and analyse pull distributions

        Parameters
        ----------
        model : str
            Model XML filename (without .xml extension)
        obs : str, optional
            Input observation definition XML filename
        trials : int, optional
            Number of trials
        duration : float, optional
            Observation duration (sec)
        ra : float, optional
            Right Ascension of pointing (deg)
        dec : float, optional
            Declination of pointing (deg)
        rad : float, optional
            Simulation radius (deg)
        emin : float, optional
            Minimum energy (TeV)
        emax : float, optional
            Maximum energy (TeV)
        enumbins : int, optional
            Number of energy bins (0 for unbinned analysis)
        edisp : bool, optional
            Use energy dispersion?
        debug : bool, optional
            Enable debugging?
        """
        # Generate pull distribution
        outfile = generate_pull_distribution(model, obs=obs, trials=trials,
                                             duration=duration,
                                             ra=ra, dec=dec, rad=rad,
                                             emin=emin, emax=emax,
                                             enumbins=enumbins,
                                             edisp=edisp, debug=debug)

        # Analyse pull distribution
        self.results = analyse_pull_distribution(outfile)

        # Return
        return

    # Test parameter result
    def test(self, name, lim_mean=0.4, lim_std=0.2):
        """
        Test one parameter

        Parameters
        ----------
        name : str
            Parameter name
        lim_mean : float, optional
            Limit for mean value
        lim_std : float, optional
            Limit for standard deviation
        """
        # Set minima and maximum
        mean_min = -lim_mean
        mean_max = +lim_mean
        std_min  = 1.0-lim_std
        std_max  = 1.0+lim_std
        
        # Test mean
        mean  = self.results[name]['mean']
        valid = (mean >= mean_min) and (mean <= mean_max)
        text  = 'Mean %.5f of %s should be within [%.2f,%.2f] range' % \
                (mean, name, mean_min, mean_max)
        self.test_assert(valid, text)

        # Test standard deviation
        std   = self.results[name]['std']
        valid = (std >= std_min) and (std <= std_max)
        text  = 'Standard deviation %.5f of %s should be within [%.2f,%.2f]' \
                ' range' % (std, name, std_min, std_max)
        self.test_assert(valid, text)

        # Return
        return

    # Test background model
    def bgd(self):
        """
        Test background model

        The standard background model is tested for an observation duration
        of 50 hours to verify the numerical accuracy of the background model
        at sufficiently good precision. Most analysis relies on the numerical
        accuracy of the background model, hence it's important to assure that
        the model is indeed accurate.
        """
        self.pull('data/sciver/bgd', duration=180000.0)
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test power law model
    def spec_plaw(self):
        """
        Test power law model
        """
        self.pull('data/sciver/crab_plaw')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test power law model with energy dispersion
    def spec_plaw_edisp(self):
        """
        Test power law model with energy dispersion
        """
        self.pull('data/sciver/crab_plaw_edisp', edisp=True)
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test power law model with stacked analysis
    def spec_plaw_stacked(self):
        """
        Test power law model with stacked analysis
        """
        self.pull('data/sciver/crab_plaw_stacked',
                  obs='data/sciver/obs_stacked.xml',
                  emin=0.020, emax=100.0, enumbins=40)
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_BackgroundModel_Prefactor')
        self.test('Pull_BackgroundModel_Index')
        return

    # Test power law 2 model
    def spec_plaw2(self):
        """
        Test power law 2 model
        """
        self.pull('data/sciver/crab_plaw2')
        self.test('Pull_Crab_PhotonFlux')
        self.test('Pull_Crab_Index')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test smoothly broken power law model
    def spec_smoothbplaw(self):
        """
        Test smoothly broken power law model
        """
        self.pull('data/sciver/crab_smoothbplaw')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index1')
        self.test('Pull_Crab_Index2')
        self.test('Pull_Crab_BreakEnergy')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test exponentially cut off power law model
    def spec_eplaw(self):
        """
        Test exponentially cut off power law model
        """
        self.pull('data/sciver/crab_eplaw')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Crab_CutoffEnergy')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test super exponentially cut off power law model
    def spec_supeplaw(self):
        """
        Test super exponentially cut off power law model
        """
        self.pull('data/sciver/crab_supeplaw')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index1')
        self.test('Pull_Crab_Index2')
        self.test('Pull_Crab_CutoffEnergy')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test log parabola model
    def spec_logparabola(self):
        """
        Test log parabola model
        """
        self.pull('data/sciver/crab_logparabola')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Crab_Curvature')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test Gaussian model
    def spec_gauss(self):
        """
        Test Gaussian model
        """
        self.pull('data/sciver/crab_gauss')
        self.test('Pull_Crab_Normalization')
        self.test('Pull_Crab_Mean')
        self.test('Pull_Crab_Sigma')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test file function model
    def spec_filefct(self):
        """
        Test file function model
        """
        self.pull('data/sciver/crab_filefct')
        self.test('Pull_Crab_Normalization')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test nodes model
    def spec_nodes(self):
        """
        Test nodes model
        """
        self.pull('data/sciver/crab_nodes')
        self.test('Pull_Crab_Intensity0')
        self.test('Pull_Crab_Intensity1')
        self.test('Pull_Crab_Intensity2')
        self.test('Pull_Crab_Intensity3')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test point source model
    def spat_ptsrc(self):
        """
        Test point source model
        """
        self.pull('data/sciver/crab_ptsrc')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Crab_RA')
        self.test('Pull_Crab_DEC')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test radial disk model
    def spat_rdisk(self):
        """
        Test radial disk model
        """
        self.pull('data/sciver/crab_rdisk')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Crab_RA')
        self.test('Pull_Crab_DEC')
        self.test('Pull_Crab_Radius')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test radial Gaussian model
    def spat_rgauss(self):
        """
        Test radial Gaussian model
        """
        self.pull('data/sciver/crab_rgauss')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Crab_RA')
        self.test('Pull_Crab_DEC')
        self.test('Pull_Crab_Sigma')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test radial shell model
    def spat_rshell(self):
        """
        Test radial shell model
        """
        self.pull('data/sciver/crab_rshell')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Crab_RA')
        self.test('Pull_Crab_DEC')
        self.test('Pull_Crab_Radius')
        self.test('Pull_Crab_Width')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test elliptical disk model
    def spat_edisk(self):
        """
        Test elliptical disk model
        """
        self.pull('data/sciver/crab_edisk')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Crab_RA')
        self.test('Pull_Crab_DEC')
        self.test('Pull_Crab_PA')
        self.test('Pull_Crab_MinorRadius')
        self.test('Pull_Crab_MajorRadius')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test elliptical Gaussian model
    def spat_egauss(self):
        """
        Test elliptical Gaussian model
        """
        self.pull('data/sciver/crab_egauss')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Crab_RA')
        self.test('Pull_Crab_DEC')
        self.test('Pull_Crab_PA')
        self.test('Pull_Crab_MinorRadius')
        self.test('Pull_Crab_MajorRadius')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test diffuse isotropic model
    def spat_const(self):
        """
        Test diffuse isotropic model
        """
        self.pull('data/sciver/crab_const')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test diffuse map model
    def spat_map(self):
        """
        Test diffuse map model
        """
        self.pull('data/sciver/crab_map', ra=201.3651, dec=-43.0191)
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test diffuse map model (small ROI)
    def spat_map_roi(self):
        """
        Test diffuse map model (small ROI)

        Note that the prefactor seems here a bit biased, which could relate
        to a possible uncertainty in the flux evaluation. This needs to be
        investigated further.
        """
        self.pull('data/sciver/crab_map_roi', ra=201.3651, dec=-43.0191, rad=1.5)
        self.test('Pull_Crab_Prefactor', lim_mean=0.45) # Accept a small bias
        self.test('Pull_Crab_Index')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test diffuse map model (not normalized and scaled)
    def spat_map_nn(self):
        """
        Test diffuse map model (not normalized and scaled)
        """
        self.pull('data/sciver/crab_map_nn', ra=201.3651, dec=-43.0191, rad=1.5)
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return

    # Test diffuse cube model
    def spat_cube(self):
        """
        Test diffuse cube model
        """
        self.pull('data/sciver/crab_cube')
        self.test('Pull_Crab_Prefactor')
        self.test('Pull_Crab_Index')
        self.test('Pull_Background_Prefactor')
        self.test('Pull_Background_Index')
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Allocate test suite container
    suites = gammalib.GTestSuites('ctools science verification')

    # Allocate test suite and append it to the container
    suite_sciver = sciver()

    # Setup test suit
    suite_sciver.set()

    # Append test suite to container
    suites.append(suite_sciver)

    # Create pfiles directory
    try:
        os.mkdir('pfiles')
    except:
        pass

    # Copy ctools parameter files into pfiles directory
    os.system('cp -r ../src/*/*.par pfiles/')

    # Set PFILES environment variable
    os.environ['PFILES'] = 'pfiles'

    # Run test suite
    success = suites.run()

    # Save test results
    suites.save('reports/sciver.xml')

    # Set return code
    if success:
        rc = 0
    else:
        rc = 1

    # Exit with return code
    sys.exit(rc)
