#! /usr/bin/env python
# ==========================================================================
# Perform unbinned in-memory analysis of simulated CTA data
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
def run_pipeline(obs, ra=83.63, dec=22.01, rad=3.0,
                 emin=0.1, emax=100.0,
                 tmin=0.0, tmax=0.0,
                 debug=False):
    """
    Simulation and binned analysis pipeline

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    ra : float, optional
        Right Ascension of Region of Interest centre (deg)
    dec : float, optional
        Declination of Region of Interest centre (deg)
    rad : float, optional
        Radius of Region of Interest (deg)
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    tmin : float, optional
        Start time (s)
    tmax : float, optional
        Stop time (s)
    debug : bool, optional
        Debug function
    """
    # Simulate events
    sim = ctools.ctobssim(obs)
    sim['debug'] = debug
    sim.run()

    # Select events
    select = ctools.ctselect(sim.obs())
    select['ra']    = ra
    select['dec']   = dec
    select['rad']   = rad
    select['emin']  = emin
    select['emax']  = emax
    select['tmin']  = tmin
    select['tmax']  = tmax
    select['debug'] = debug
    select.run()

    # Perform maximum likelihood fitting
    like = ctools.ctlike(select.obs())
    like['debug'] = True # Switch this always on for results in console
    like.run()

    # Return
    return


# =============================== #
# Run unbinned in-memory pipeline #
# =============================== #
def pipeline_unbinned_mem():
    """
    Run unbinned in-memory pipeline
    """
    # Set usage string
    usage = 'pipeline_unbinned_mem.py [-d datadir]'

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
    run_pipeline(obs)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Run unbinned in-memory pipeline
    pipeline_unbinned_mem()
