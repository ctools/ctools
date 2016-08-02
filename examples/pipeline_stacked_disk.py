#! /usr/bin/env python
# ==========================================================================
# Perform stacked analysis of simulated CTA data.
#
# Copyright (C) 2015-2016 Juergen Knoedlseder
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
import cscripts


# ================================ #
# Simulation and analysis pipeline #
# ================================ #
def run_pipeline(obs, ra=83.63, dec=22.01, emin=0.1, emax=100.0,
                 enumbins=20, nxpix=200, nypix=200, binsz=0.02,
                 coordsys='CEL', proj='CAR',
                 model='data/crab.xml',
                 caldb='prod2', irf='South_0.5h',
                 debug=False):
    """
    Simulation and stacked analysis pipeline

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    ra : float, optional
        Right Ascension of counts cube centre (deg)
    dec : float, optional
        Declination of Region of counts cube centre (deg)
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    enumbins : int, optional
        Number of energy bins
    nxpix : int, optional
        Number of pixels in X axis
    nypix : int, optional
        Number of pixels in Y axis
    binsz : float, optional
        Pixel size (deg)
    coordsys : str, optional
        Coordinate system
    proj : str, optional
        Coordinate projection
    model : str, optional
        Model definition XML file
    caldb : str, optional
        Calibration database path
    irf : str, optional
        Instrument response function
    debug : bool, optional
        Debug function
    """
    # Simulate events
    sim = ctools.ctobssim(obs)
    sim['debug']     = debug
    sim['outevents'] = 'obs.xml'
    sim.execute()

    # Bin events into counts map
    bin = ctools.ctbin()
    bin['inobs']    = 'obs.xml'
    bin['outcube']  = 'cntcube.fits'
    bin['ebinalg']  = 'LOG'
    bin['emin']     = emin
    bin['emax']     = emax
    bin['enumbins'] = enumbins
    bin['nxpix']    = nxpix
    bin['nypix']    = nypix
    bin['binsz']    = binsz
    bin['coordsys'] = coordsys
    bin['proj']     = proj
    bin['xref']     = ra
    bin['yref']     = dec
    bin['debug']    = debug
    bin.execute()

    # Create exposure cube
    expcube = ctools.ctexpcube()
    expcube['inobs']    = 'obs.xml'
    expcube['incube']   = 'cntcube.fits'
    expcube['outcube']  = 'expcube.fits'
    expcube['caldb']    = caldb
    expcube['irf']      = irf
    expcube['ebinalg']  = 'LOG'
    expcube['emin']     = emin
    expcube['emax']     = emax
    expcube['enumbins'] = enumbins
    expcube['nxpix']    = nxpix
    expcube['nypix']    = nypix
    expcube['binsz']    = binsz
    expcube['coordsys'] = coordsys
    expcube['proj']     = proj
    expcube['xref']     = ra
    expcube['yref']     = dec
    expcube['debug']    = debug
    expcube.execute()

    # Create PSF cube
    psfcube = ctools.ctpsfcube()
    psfcube['inobs']    = 'obs.xml'
    psfcube['incube']   = 'NONE'
    psfcube['outcube']  = 'psfcube.fits'
    psfcube['caldb']    = caldb
    psfcube['irf']      = irf
    psfcube['ebinalg']  = 'LOG'
    psfcube['emin']     = emin
    psfcube['emax']     = emax
    psfcube['enumbins'] = enumbins
    psfcube['nxpix']    = 10
    psfcube['nypix']    = 10
    psfcube['binsz']    = 1.0
    psfcube['coordsys'] = coordsys
    psfcube['proj']     = proj
    psfcube['xref']     = ra
    psfcube['yref']     = dec
    psfcube['debug']    = debug
    psfcube.execute()

    # Create background cube
    bkgcube = ctools.ctbkgcube()
    bkgcube['inobs']    = 'obs.xml'
    bkgcube['inmodel']  = model
    bkgcube['incube']   = 'NONE'
    bkgcube['outcube']  = 'bkgcube.fits'
    bkgcube['outmodel'] = 'model_bkg.xml'
    bkgcube['caldb']    = caldb
    bkgcube['irf']      = irf
    bkgcube['ebinalg']  = 'LOG'
    bkgcube['emin']     = emin
    bkgcube['emax']     = emax
    bkgcube['enumbins'] = enumbins
    bkgcube['nxpix']    = 10
    bkgcube['nypix']    = 10
    bkgcube['binsz']    = 1.0
    bkgcube['coordsys'] = coordsys
    bkgcube['proj']     = proj
    bkgcube['xref']     = ra
    bkgcube['yref']     = dec
    bkgcube['debug']    = debug
    bkgcube.execute()

    # Perform maximum likelihood fitting
    like = ctools.ctlike()
    like['inobs']    = 'cntcube.fits'
    like['inmodel']  = 'model_bkg.xml'
    like['outmodel'] = 'fit_results.xml'
    like['expcube']  = 'expcube.fits'
    like['psfcube']  = 'psfcube.fits'
    like['bkgcube']  = 'bkgcube.fits'
    like['caldb']    = caldb
    like['irf']      = irf
    like['debug']    = True # Switch this always on for results in console
    like.execute()

    # Return
    return


# ====================================================== #
# Run stacked pipeline with intermediate results on disk #
# ====================================================== #
def pipeline_stacked_disk():
    """
    Run stacked pipeline with intermediate results on disk
    """
    # Set usage string
    usage = 'pipeline_stacked_disk.py [-d datadir]'

    # Set default options
    options = [{'option': '-d', 'value': 'data'}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    datadir = options[0]['value']

    # Setup observations
    obs = cscripts.obsutils.set_observations(83.63, 22.01, 5.0, 0.0, 180.0,
                                             0.1, 100.0, 'South_0.5h', 'prod2',
                                             pattern='four', offset=1.5)

    # Setup model
    obs.models(gammalib.GModels(datadir+'/crab.xml'))

    # Run analysis pipeline
    run_pipeline(obs, model=datadir+'/crab.xml', enumbins=10, nxpix=40,
                 nypix=40, binsz=0.1)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Run stacked pipeline with intermediate results on disk
    pipeline_stacked_disk()
