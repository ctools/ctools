# ==========================================================================
# This script provides a number of functions that are useful for handling
# CTA observations.
#
# Copyright (C) 2011-2015 Juergen Knoedlseder
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
import math


# ===================== #
# Simulate observations #
# ===================== #
def sim(obs, log=False, debug=False, chatter=2, edisp=False, seed=0, nbins=0,
        binsz=0.05, npix=200, proj="TAN", coord="GAL"):
    """
    Simulate events for all observations in the container.

    Parameters:
     obs   - Observation container
    Keywords:
     log   - Create log file(s)
     debug - Create console dump?
     edisp - Apply energy dispersion?
     seed  - Seed value for simulations (default: 0)
     nbins - Number of energy bins (default: 0=unbinned)
     binsz - Pixel size for binned simulation (deg/pixel)
     npix  - Number of pixels in X and Y for binned simulation
    """

    # Allocate ctobssim application and set parameters
    sim = ctools.ctobssim(obs)
    sim["seed"] = seed
    sim["edisp"] = edisp

    # Optionally open the log file
    if log:
        sim.logFileOpen()

    # Optionally switch-on debugging model
    if debug:
        sim["debug"] = True

    # Set chatter level
    sim["chatter"] = chatter

    # Run ctobssim application. This will loop over all observations in the
    # container and simulation the events for each observation. Note that
    # events are not added together, they still apply to each observation
    # separately.
    sim.run()

    # Binned option?
    if nbins > 0:

        # Determine common energy boundaries for observations
        emin = None
        emax = None
        for run in sim.obs():
            run_emin = run.events().ebounds().emin().TeV()
            run_emax = run.events().ebounds().emax().TeV()
            if emin == None:
                emin = run_emin
            elif run_emin > emin:
                emin = run_emin
            if emax == None:
                emax = run_emax
            elif run_emax > emax:
                emax = run_emax

        # Allocate ctbin application and set parameters
        bin = ctools.ctbin(sim.obs())
        bin["ebinalg"] = "LOG"
        bin["emin"] = emin
        bin["emax"] = emax
        bin["enumbins"] = nbins
        bin["usepnt"] = True # Use pointing for map centre
        bin["nxpix"] = npix
        bin["nypix"] = npix
        bin["binsz"] = binsz
        bin["coordsys"] = coord
        bin["proj"] = proj

        # Optionally open the log file
        if log:
            bin.logFileOpen()

        # Optionally switch-on debugging model
        if debug:
            bin["debug"] = True

        # Set chatter level
        bin["chatter"] = chatter

        # Run ctbin application. This will loop over all observations in
        # the container and bin the events in counts maps
        bin.run()

        # Make a deep copy of the observation that will be returned
        # (the ctbin object will go out of scope one the function is
        # left)
        obs = bin.obs().copy()

    else:

        # Make a deep copy of the observation that will be returned
        # (the ctobssim object will go out of scope one the function is
        # left)
        obs = sim.obs().copy()

    # Delete the simulation
    del sim

    # Return observation container
    return obs


# ================ #
# Fit observations #
# ================ #
def fit(obs, log=False, debug=False, chatter=2, edisp=False):
    """
    Perform maximum likelihood fitting of observations in the container.

    Parameters:
     obs   - Observation container
    Keywords:
     log     - Create log file(s)
     debug   - Create screen dump
     chatter - Chatter level
     edisp   - Apply energy dispersion?
    """
    # Allocate ctlike application
    like = ctools.ctlike(obs)

    # Optionally open the log file
    if log:
        like.logFileOpen()

    # Optionally switch-on debugging model
    if debug:
        like["debug"] = True

    # Set chatter level
    like["chatter"] = chatter

    # Optionally apply energy dispersion
    like["edisp"] = edisp

    # Run ctlike application.
    like.run()

    # Return observations
    return like


# ============================================================== #
# Fit observations and determine errors using likelihood profile #
# ============================================================== #
def cterror(obs, srcname, log=False, debug=False, chatter=2):
    """
    Perform maximum likelihood fitting of observations in the container.

    Parameters:
     obs     - Observation container
     srcname - Source name
    Keywords:
     log     - Create log file(s)
     debug   - Create screen dump
     chatter - Chatter level
     edisp   - Apply energy dispersion?
    """
    # Allocate cterror application
    error = ctools.cterror(obs)

    # Set cterror parameters
    error["srcname"] = srcname

    # Optionally open the log file
    if log:
        error.logFileOpen()

    # Optionally switch-on debugging model
    if debug:
        error["debug"] = True

    # Set chatter level
    error["chatter"] = chatter

    # Run cterror application.
    error.run()

    # Return observations
    return error


# ================= #
# Create counts map #
# ================= #
def cntmap(obs, proj="TAN", coord="GAL", xval=0.0, yval=0.0, \
           binsz=0.05, nxpix=200, nypix=200, \
           outname="cntmap.fits"):
    """
    Creates a counts map by combining the events of all observations.
    The counts map will be a summed map over all energies.

    Parameters:
     obs     - Observation container
    Keywords:
     proj    - Projection type (e.g. TAN, CAR, STG, ...) (default: TAN)
     coord   - Coordinate type (GAL, CEL) (default: GAL)
     xval    - Reference longitude value [deg] (default: 0.0)
     yval    - Reference latitude value [deg] (default: 0.0)
     binsz   - Pixel size [deg/pixel] (default: 0.05)
     nxpix   - Number of pixels in X direction (default: 200)
     nypix   - Number of pixels in Y direction (default: 200)
     outname - Counts map FITS filename (default: cntmap.fits)
    """
    # Allocate counts map
    map = gammalib.GSkyMap(proj, coord, xval, yval, -binsz, binsz, nxpix, nypix, 1)

    # Set maximum pixel number
    maxpixel = nxpix * nypix

    # Fill all observations
    for run in obs:

        # Loop over all events
        for event in run.events():

            # Determine sky pixel
            skydir = event.dir().dir()
            pixel  = map.dir2inx(skydir)

            # Set pixel
            if pixel < maxpixel:
                map[pixel] += event.counts()

    # Save sky map. The clobber flag is set to True, so any existing FITS
    # file will be overwritten.
    map.save(outname, True)

    # Return counts map
    return map


# ================ #
# Create model map #
# ================ #
def modmap(obs, eref=0.1, proj="TAN", coord="GAL", xval=0.0, yval=0.0, \
           binsz=0.05, nxpix=200, nypix=200, \
           outname="modmap.fits"):
    """
    Make model map for a given reference energy by combining all observations.
    The model map will be evaluated for a given reference energy 'eref' and will
    be given in units of [counts/(sr MeV s)].

    Parameters:
     obs     - Observation container
    Keywords:
     eref    - Reference energy for which model is created [TeV] (default: 0.1)
     proj    - Projection type (e.g. TAN, CAR, STG, ...) (default: TAN)
     coord   - Coordinate type (GAL, CEL) (default: GAL)
     xval    - Reference longitude value [deg] (default: 0.0)
     yval    - Reference latitude value [deg] (default: 0.0)
     binsz   - Pixel size [deg/pixel] (default: 0.05)
     nxpix   - Number of pixels in X direction (default: 200)
     nypix   - Number of pixels in Y direction (default: 200)
     outname - Model map FITS filename (default: modmap.fits)
    """
    # Allocate model map
    map = gammalib.GSkyMap(proj, coord, xval, yval, -binsz, binsz, nxpix, nypix, 1)

    # Set reference energy, time and direction. The time is not initialised and is
    # in fact not used (as the IRF is assumed to be time independent for now).
    # The sky direction is set later using the pixel values.
    energy  = gammalib.GEnergy()
    time    = gammalib.GTime()
    instdir = gammalib.GCTAInstDir()
    energy.TeV(eref)

    # Loop over all map pixels
    for pixel in range(map.npix()):

        # Get sky direction
        skydir = map.inx2dir(pixel)
        instdir.dir(skydir)

        # Create event atom for map pixel
        atom = gammalib.GCTAEventAtom()
        atom.dir(instdir)
        atom.energy(energy)
        atom.time(time)

        # Initialise model value
        value = 0.0

        # Loop over all observations
        for run in obs:
            value += obs.models().eval(atom, run)

        # Set map value
        map[pixel] = value

    # Save sky map
    map.save(outname, True)

    # Return model map
    return map


# ===================================== #
# Set one CTA observation (old version) #
# ===================================== #
def set(pntdir, tstart=0.0, duration=1800.0, deadc=0.95, \
        emin=0.1, emax=100.0, rad=5.0, \
        irf="South_50h", caldb="prod2"):
    """
    Obsolete function, use set_obs instead.
    """
    # Print warning
    print("Warning: obsutils.set is obsolete, use obsutils.set_obs instead.")

    # Call new function
    obs_cta = set_obs(pntdir, tstart=tstart, duration=duration, deadc=deadc, \
                      emin=emin, emax=emax, rad=rad, \
                      irf=irf, caldb=caldb)

    # Return CTA observation
    return obs_cta


# ======================= #
# Set one CTA observation #
# ======================= #
def set_obs(pntdir, tstart=0.0, duration=1800.0, deadc=0.95, \
            emin=0.1, emax=100.0, rad=5.0, \
            irf="South_50h", caldb="prod2", id="000000", instrument="CTA"):
    """
    Returns a single CTA observation containing an empty CTA event list.
    By looping over this function you can add CTA observations to the
    observation container.

    Parameters:
     pntdir   - Pointing direction [GSkyDir]
    Keywords:
     tstart     - Start time [seconds] (default: 0.0)
     duration   - Duration of observation [seconds] (default: 1800.0)
     deadc      - Deadtime correction factor (default: 0.95)
     emin       - Minimum event energy [TeV] (default: 0.1)
     emax       - Maximum event energy [TeV] (default: 100.0)
     rad        - ROI radius used for analysis [deg] (default: 5.0)
     irf        - Instrument response function (default: cta_dummy_irf)
     caldb      - Calibration database path (default: "dummy")
     id         - Run identifier (default: "000000")
     instrument - Intrument (default: "CTA")
    """
    # Allocate CTA observation
    obs_cta = gammalib.GCTAObservation(instrument)

    # Set calibration database
    db = gammalib.GCaldb()
    if (gammalib.dir_exists(caldb)):
        db.rootdir(caldb)
    else:
        db.open("cta", caldb)

    # Set pointing direction
    pnt = gammalib.GCTAPointing()
    pnt.dir(pntdir)
    obs_cta.pointing(pnt)

    # Set ROI
    roi     = gammalib.GCTARoi()
    instdir = gammalib.GCTAInstDir()
    instdir.dir(pntdir)
    roi.centre(instdir)
    roi.radius(rad)

    # Set GTI
    gti = gammalib.GGti()
    gti.append(gammalib.GTime(tstart), gammalib.GTime(tstart+duration))

    # Set energy boundaries
    ebounds = gammalib.GEbounds(gammalib.GEnergy(emin, "TeV"), \
                                gammalib.GEnergy(emax, "TeV"))

    # Allocate event list
    events = gammalib.GCTAEventList()
    events.roi(roi)
    events.gti(gti)
    events.ebounds(ebounds)
    obs_cta.events(events)

    # Set instrument response
    obs_cta.response(irf, db)

    # Set ontime, livetime, and deadtime correction factor
    obs_cta.ontime(duration)
    obs_cta.livetime(duration*deadc)
    obs_cta.deadc(deadc)
    obs_cta.id(id)

    # Return CTA observation
    return obs_cta


# ============================ #
# Set list of CTA observations #
# ============================ #
def set_obs_list(obsdeflist, tstart=0.0, duration=1800.0, deadc=0.95, \
        emin=0.1, emax=100.0, rad=5.0, \
        irf="South_50h", caldb="prod2"):
    """
    Returns an observation container filled with a list of CTA observations.
    The list is defined by the obsdeflist parameter which is a dictionnary
    containing the mandatory keywords 'ra' and 'dec' that specify the
    pointing direction for a given observation. Optional keyword give control
    over other observation proposerties, such as duration, deadtime correction,
    energy range, etc. If an optional keyword is not specified, the function
    keyword is used instead.

    Parameters:
     obsdeflist - Observation definition list [{'ra': x.xx, 'dec': x.xx}]
                  The directory can take the following optional keywords:
                  - duration: Duration of CTA observation [seconds]
                  - deadc: Deadtime correction factor
                  - emin: Minimum event energy [TeV]
                  - emax: Maximum event energy [TeV]
                  - rad: ROI radius used for analysis [deg]
                  - irf: Instrument response function
                  - caldb: Calibration database path
                  Optional keywords overwrite keywords specified in the
                  function call.
    Keywords:
     tstart     - Start time [seconds] (default: 0.0)
     duration   - Duration of one CTA observation [seconds] (default: 1800.0)
     deadc      - Deadtime correction factor (default: 0.95)
     emin       - Minimum event energy [TeV] (default: 0.1)
     emax       - Maximum event energy [TeV] (default: 100.0)
     rad        - ROI radius used for analysis [deg] (default: 5.0)
     irf        - Instrument response function (default: cta_dummy_irf)
     caldb      - Calibration database path (default: "dummy")
    """
    # Initialise empty observation container
    obs = gammalib.GObservations()

    # Initialise first time and identifier
    obs_start = tstart
    obs_id    = 1

    # Loop over observation definition list
    for obsdef in obsdeflist:

        # Set pointing direction
        pntdir = gammalib.GSkyDir()
        pntdir.radec_deg(obsdef['ra'], obsdef['dec'])

        # Set duration (use default if not found in definition list)
        if 'duration' in obsdef:
            obs_duration = obsdef['duration']
        else:
            obs_duration = duration

        # Set emin (use default if not found in definition list)
        if 'emin' in obsdef:
            obs_emin = obsdef['emin']
        else:
            obs_emin = emin

        # Set emax (use default if not found in definition list)
        if 'emax' in obsdef:
            obs_emax = obsdef['emax']
        else:
            obs_emax = emax

        # Set radius (use default if not found in definition list)
        if 'rad' in obsdef:
            obs_rad = obsdef['rad']
        else:
            obs_rad = rad

        # Set deadc (use default if not found in definition list)
        if 'deadc' in obsdef:
            obs_deadc = obsdef['deadc']
        else:
            obs_deadc = deadc

        # Set caldb (use default if not found in definition list)
        if 'caldb' in obsdef:
            obs_caldb = obsdef['caldb']
        else:
            obs_caldb = caldb

        # Set irf (use default if not found in definition list)
        if 'irf' in obsdef:
            obs_irf = obsdef['irf']
        else:
            obs_irf = irf

        # Generate identifier string
        id = "%6.6d" % obs_id

        # Set CTA observation
        obs_cta = set_obs(pntdir, tstart=obs_start, duration=obs_duration, \
                          deadc=obs_deadc, emin=obs_emin, emax=obs_emax, \
                          rad=obs_rad, irf=obs_irf, caldb=obs_caldb, id=id)

        # Append to container
        obs.append(obs_cta)

        # Update start time and identifier
        obs_start += obs_duration
        obs_id    += 1

    # Return observation container
    return obs


# ============================ #
# Set CTA observation patterns #
# ============================ #
def set_obs_patterns(pattern, ra=83.6331, dec=22.0145, offset=1.5):
    """
    Sets a number of standard patterns.

    Parameters:
     pattern - Observation pattern. Possible options are:
               - "single": single pointing
               - "four": four pointings 'offset' around pattern centre

    Keywords:
     ra      - Right Ascension of pattern centre [deg] (default: 83.6331)
     dec     - Declination of pattern centre [deg] (default: 22.0145)
     offset  - Offset from pattern centre [deg] (default: 1.5)
    """
    # Initialise observation definition list
    obsdeflist = []

    # Add patterns
    if pattern == "single":
        obsdef = {'ra': ra, 'dec': dec}
        obsdeflist.append(obsdef)
    elif pattern == "four":
        # Set pattern centre
        centre = gammalib.GSkyDir()
        centre.radec_deg(ra, dec)

        # Append pointings
        for phi in [0.0, 90.0, 180.0, 270.0]:
            pntdir = centre.copy()
            pntdir.rotate_deg(phi, offset)
            obsdeflist.append({'ra': pntdir.ra_deg(), 'dec': pntdir.dec_deg()})
    else:
        print("Warning: Observation pattern '"+str(pattern)+"' not recognized.")

    # Return observation definition list
    return obsdeflist
