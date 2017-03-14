# ==========================================================================
# Utility functions for observation handling
#
# Copyright (C) 2011-2017 Juergen Knoedlseder
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
import math
import gammalib
import ctools


# ===================== #
# Simulate observations #
# ===================== #
def sim(obs, log=False, debug=False, chatter=2, edisp=False, seed=0,
        emin=None, emax=None, nbins=0, addbounds=False,
        binsz=0.05, npix=200, proj='TAN', coord='GAL'):
    """
    Simulate events for all observations in the container

    Simulate events for all observations using ctobssim. If the number of
    energy bins is positive, the events are binned into a counts cube using
    ctbin. If multiple observations are simulated, the counts cube is a
    stacked cube and the corresponding response cubes are computed using
    ctexpcube, ctpsfcube, ctbkgcube and optionally ctedispcube. The response
    cubes are attached to the first observation in the container, which
    normally is the observation with the counts cube.

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container without events
    log : bool, optional
        Create log file(s)
    debug : bool, optional
        Create console dump?
    chatter : int, optional
        Chatter level
    edisp : bool, optional
        Apply energy dispersion?
    seed : int, optional
        Seed value for simulations
    emin : float, optional
        Minimum energy of counts cube for binned (TeV)
    emax : float, optional
        Maximum energy of counts cube for binned (TeV)
    nbins : int, optional
        Number of energy bins (0=unbinned)
    addbounds : bool, optional
        Add boundaries at observation energies
    binsz : float, optional
        Pixel size for binned simulation (deg/pixel)
    npix : int, optional
        Number of pixels in X and Y for binned simulation
    proj : str, optional
        Projection for binned simulation
    coord : str, optional
        Coordinate system for binned simulation

    Returns
    -------
    obs : `~gammalib.GObservations`
        Observation container filled with simulated events
    """
    # Allocate ctobssim application and set parameters
    sim = ctools.ctobssim(obs)
    sim['seed']    = seed
    sim['edisp']   = edisp
    sim['chatter'] = chatter
    sim['debug']   = debug

    # Optionally open the log file
    if log:
        sim.logFileOpen()

    # Run ctobssim application. This will loop over all observations in the
    # container and simulation the events for each observation. Note that
    # events are not added together, they still apply to each observation
    # separately.
    sim.run()

    # Binned option?
    if nbins > 0:

        # If energy boundaries are not given then determine the minimum and
        # the maximum energies from all observations and use these values
        # as energy boundaries. The energy boundaries are given in TeV.
        if emin == None or emax == None:
            emin = 1.0e30
            emax = 0.0
            for run in sim.obs():
                emin = min(run.events().ebounds().emin().TeV(), emin)
                emax = max(run.events().ebounds().emax().TeV(), emax)

        # Allocate ctbin application and set parameters
        bin = ctools.ctbin(sim.obs())
        bin['ebinalg']  = 'LOG'
        bin['emin']     = emin
        bin['emax']     = emax
        bin['enumbins'] = nbins
        bin['usepnt']   = True # Use pointing for map centre
        bin['nxpix']    = npix
        bin['nypix']    = npix
        bin['binsz']    = binsz
        bin['coordsys'] = coord
        bin['proj']     = proj
        bin['chatter']  = chatter
        bin['debug']    = debug

        # Optionally open the log file
        if log:
            bin.logFileOpen()

        # Run ctbin application. This will loop over all observations in
        # the container and bin the events in counts maps
        bin.run()

        # If we have multiple input observations then create stacked response
        # cubes and append them to the observation
        if len(sim.obs()) > 1:

            # Get stacked response (use pointing for map centre)
            response = get_stacked_response(sim.obs(), None, None,
                                            binsz=binsz, nxpix=npix, nypix=npix,
                                            emin=emin, emax=emax, enumbins=nbins,
                                            edisp=edisp,
                                            coordsys=coord, proj=proj,
                                            addbounds=addbounds,
                                            log=log, debug=debug,
                                            chatter=chatter)

            # Set stacked response
            if edisp:
                bin.obs()[0].response(response['expcube'],
                                      response['psfcube'],
                                      response['edispcube'],
                                      response['bkgcube'])
            else:
                bin.obs()[0].response(response['expcube'],
                                      response['psfcube'],
                                      response['bkgcube'])

            # Set new models
            bin.obs().models(response['models'])

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


# ======================= #
# Set one CTA observation #
# ======================= #
def set_obs(pntdir, tstart=0.0, duration=1800.0, deadc=0.98, \
            emin=0.1, emax=100.0, rad=5.0, \
            irf='South_50h', caldb='prod2', obsid='000000'):
    """
    Set a single CTA observation
    
    The function sets a single CTA observation containing an empty CTA
    event list. By looping over this function CTA observations can be
    added to the observation container.

    Parameters
    ----------
    pntdir : `~gammalib.GSkyDir`
        Pointing direction
    tstart : float, optional
        Start time (s)
    duration : float, optional
        Duration of observation (s)
    deadc : float, optional
        Deadtime correction factor
    emin : float, optional
        Minimum event energy (TeV)
    emax : float, optional
        Maximum event energy (TeV)
    rad : float, optional
        ROI radius used for analysis (deg)
    irf : str, optional
        Instrument response function
    caldb : str, optional
        Calibration database path
    obsid : str, optional
        Observation identifier

    Returns
    -------
    obs : `~gammalib.GCTAObservation`
        CTA observation
    """
    # Allocate CTA observation
    obs = gammalib.GCTAObservation()

    # Set CTA calibration database
    db = gammalib.GCaldb()
    if (gammalib.dir_exists(caldb)):
        db.rootdir(caldb)
    else:
        db.open('cta', caldb)

    # Set pointing direction for CTA observation
    pnt = gammalib.GCTAPointing()
    pnt.dir(pntdir)
    obs.pointing(pnt)

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
    ebounds = gammalib.GEbounds(gammalib.GEnergy(emin, 'TeV'),
                                gammalib.GEnergy(emax, 'TeV'))

    # Allocate event list
    events = gammalib.GCTAEventList()

    # Set ROI, GTI and energy boundaries for event list
    events.roi(roi)
    events.gti(gti)
    events.ebounds(ebounds)

    # Set the event list as the events for CTA observation
    obs.events(events)

    # Set instrument response for CTA observation
    obs.response(irf, db)

    # Set ontime, livetime, and deadtime correction factor for CTA observation
    obs.ontime(duration)
    obs.livetime(duration*deadc)
    obs.deadc(deadc)
    obs.id(obsid)

    # Return CTA observation
    return obs


# ============================ #
# Set list of CTA observations #
# ============================ #
def set_obs_list(obsdeflist, tstart=0.0, duration=1800.0, deadc=0.98, \
                 emin=0.1, emax=100.0, rad=5.0, \
                 irf='South_50h', caldb='prod2'):
    """
    Returns an observation container filled with a list of CTA observations

    The list is defined by the obsdeflist parameter which is a dictionnary
    containing the mandatory keywords 'ra' and 'dec' that specify the
    pointing direction for a given observation. Optional keyword give control
    over other observation proposerties, such as duration, deadtime correction,
    energy range, etc. If an optional keyword is not specified, the function
    keyword is used instead.

    Parameters
    ----------
    obsdeflist : list of dict
        Observation definition list
    tstart : float, optional
        Start time (s)
    duration : float, optional
        Duration of observation (s)
    deadc : float, optional
        Deadtime correction factor
    emin : float, optional
        Minimum event energy (TeV)
    emax : float, optional
        Maximum event energy (TeV)
    rad : float, optional
        ROI radius used for analysis (deg)
    irf : str, optional
        Instrument response function
    caldb : str, optional
        Calibration database path

    Returns
    -------
    obs : `~gammalib.GObservations`
        Observation container filled with CTA observation
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

        # Generate identifier string
        obsid = '%6.6d' % obs_id

        # Set one CTA observation
        obs_cta = set_obs(pntdir,
                          tstart=obs_start,
                          duration=obsdef.setdefault('duration', duration),
                          deadc=obsdef.setdefault('deadc', deadc),
                          emin=obsdef.setdefault('emin', emin),
                          emax=obsdef.setdefault('emax', emax),
                          rad=obsdef.setdefault('rad', rad),
                          irf=obsdef.setdefault('irf', irf),
                          caldb=obsdef.setdefault('caldb', caldb),
                          obsid=obsid)

        # Append to container
        obs.append(obs_cta)

        # Update start time and identifier
        obs_start += duration
        obs_id    += 1

    # Return observation container
    return obs


# ============================ #
# Set CTA observation patterns #
# ============================ #
def set_obs_patterns(pattern, ra=83.6331, dec=22.0145, offset=1.5):
    """
    Sets a number of standard patterns
    
    Parameters
    ----------
    pattern : str
        Observation pattern ("single", "four")
    ra : float, optional
        Right Ascension of pattern centre (deg)
    dec : float, optional
        Declination of pattern centre (deg)
    offset : float, optional
        Offset from pattern centre (deg)

    Returns
    -------
    obsdeflist : list
        Observation definition list
    """
    # Initialise observation definition list
    obsdeflist = []

    # If the pattern is a single observation then append the Right Ascension
    # and Declination to the observation definition list
    if pattern == 'single':
        obsdef = {'ra': ra, 'dec': dec}
        obsdeflist.append(obsdef)

    # ... otherwise, if the pattern is four observations then append four
    # observations offset by a certain amount from the pattern centre to the
    # observation definition list
    elif pattern == 'four':

        # Set pattern centre
        centre = gammalib.GSkyDir()
        centre.radec_deg(ra, dec)

        # Append pointings
        for phi in [0.0, 90.0, 180.0, 270.0]:
            pntdir = centre.copy()
            pntdir.rotate_deg(phi, offset)
            obsdeflist.append({'ra': pntdir.ra_deg(), 'dec': pntdir.dec_deg()})

    # ... otherwise we have an unknown pattern
    else:
        print('Warning: Observation pattern "'+str(pattern)+'" not recognized.')

    # Return observation definition list
    return obsdeflist


# ====================================================== #
# Set observation container filled with CTA observations #
# ====================================================== #
def set_observations(ra, dec, rad, tstart, duration, emin, emax, irf, caldb,
                     deadc=0.98, pattern='single', offset=1.5):
    """
    Set an observation container filled with CTA observations
    
    Parameters
    ----------
    ra : float
        Right Ascension of pattern centre (deg)
    dec : float
        Declination of pattern centre (deg)
    rad : float
        ROI radius used for analysis (deg)
    tstart : float
        Start time of observation (s)
    duration : float
        Duration of each observation (s)
    emin : float, optional
        Minimum event energy (TeV)
    emax : float, optional
        Maximum event energy (TeV)
    irf : str
        Instrument response function
    caldb : str
        Calibration database path
    deadc : float, optional
        Deadtime correction factor
    pattern : str, optional
        Observation pattern ("single", "four")
    offset : float, optional
        Offset from pattern centre (deg)

    Returns
    -------
    obs : `~gammalib.GObservations()
        Observation container
    """
    # Setup observation definition list
    obsdeflist = set_obs_patterns(pattern, offset=offset, ra=ra, dec=dec)

    # Create list of observations
    obs = set_obs_list(obsdeflist, tstart=tstart, duration=duration,
                       emin=emin, emax=emax, rad=rad, irf=irf, caldb=caldb,
                       deadc=deadc)

    # Return observation container
    return obs


# ==================== #
# Get stacked response #
# ==================== #
def get_stacked_response(obs, xref, yref, binsz=0.05, nxpix=200, nypix=200,
                         emin=0.1, emax=100.0, enumbins=20, edisp=False,
                         coordsys='GAL', proj='TAN', addbounds=False,
                         log=False, debug=False, chatter=2):
    """
    Get stacked response cubes

    The number of energies bins are set to at least 30 bins per decade, and
    the "enumbins" parameter is only used if the number of bins is larger
    than 30 bins per decade.

    If the xref or yref arguments are "None" the response cube centre will be
    determined from the pointing information in the observation container.

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    xref : float
        Right Ascension or Galactic longitude of response centre (deg)
    yref : float
        Declination or Galactic latitude of response centre (deg)
    binsz : float, optional
        Pixel size (deg/pixel)
    nxpix : int, optional
        Number of pixels in X direction
    nypix : int, optional
        Number of pixels in Y direction
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    enumbins : int, optional
        Number of energy bins
    edisp : bool, optional
        Apply energy dispersion?
    coordsys : str, optional
        Coordinate system
    proj : str, optional
        Projection
    addbounds : bool, optional
        Add boundaries at observation energies
    log : bool, optional
        Create log file(s)
    debug : bool, optional
        Create console dump?
    chatter : int, optional
        Chatter level

    Returns
    -------
    result : dict
        Dictionary of response cubes
    """
    # If no xref and yref arguments have been specified then use the pointing
    # information
    if xref == None or yref == None:
        usepnt = True
    else:
        usepnt = False

    # Set number of energy bins to at least 30 per energy decade
    _enumbins = int((math.log10(emax) - math.log10(emin)) * 30.0)
    if enumbins > _enumbins:
        _enumbins = enumbins

    # Compute spatial binning for point spread function and energy dispersion
    # cubes. The spatial binning is 10 times coarser than the spatial binning
    # of the exposure and background cubes. At least 2 spatial are required.
    psf_binsz = 10.0 * binsz
    psf_nxpix = max(nxpix // 10, 2)  # Make sure result is int
    psf_nypix = max(nypix // 10, 2)  # Make sure result is int

    # Create exposure cube
    expcube = ctools.ctexpcube(obs)
    expcube['incube']    = 'NONE'
    expcube['usepnt']    = usepnt
    expcube['ebinalg']   = 'LOG'
    expcube['binsz']     = binsz
    expcube['nxpix']     = nxpix
    expcube['nypix']     = nypix
    expcube['enumbins']  = _enumbins
    expcube['emin']      = emin
    expcube['emax']      = emax
    expcube['coordsys']  = coordsys
    expcube['proj']      = proj
    expcube['addbounds'] = addbounds
    expcube['debug']     = debug
    expcube['chatter']   = chatter
    if not usepnt:
        expcube['xref'] = xref
        expcube['yref'] = yref
    if log:
        expcube.logFileOpen()
    expcube.run()

    # Create point spread function cube
    psfcube = ctools.ctpsfcube(obs)
    psfcube['incube']    = 'NONE'
    psfcube['usepnt']    = usepnt
    psfcube['ebinalg']   = 'LOG'
    psfcube['binsz']     = psf_binsz
    psfcube['nxpix']     = psf_nxpix
    psfcube['nypix']     = psf_nypix
    psfcube['enumbins']  = _enumbins
    psfcube['emin']      = emin
    psfcube['emax']      = emax
    psfcube['coordsys']  = coordsys
    psfcube['proj']      = proj
    psfcube['addbounds'] = addbounds
    psfcube['debug']     = debug
    psfcube['chatter']   = chatter
    if not usepnt:
        psfcube['xref'] = xref
        psfcube['yref'] = yref
    if log:
        psfcube.logFileOpen()
    psfcube.run()

    # Create background cube
    bkgcube = ctools.ctbkgcube(obs)
    bkgcube['incube']    = 'NONE'
    bkgcube['usepnt']    = usepnt
    bkgcube['ebinalg']   = 'LOG'
    bkgcube['binsz']     = binsz
    bkgcube['nxpix']     = nxpix
    bkgcube['nypix']     = nypix
    bkgcube['enumbins']  = _enumbins
    bkgcube['emin']      = emin
    bkgcube['emax']      = emax
    bkgcube['coordsys']  = coordsys
    bkgcube['proj']      = proj
    bkgcube['addbounds'] = addbounds
    bkgcube['debug']     = debug
    bkgcube['chatter']   = chatter
    if not usepnt:
        bkgcube['xref'] = xref
        bkgcube['yref'] = yref
    if log:
        bkgcube.logFileOpen()
    bkgcube.run()

    # If energy dispersion is requested then create energy dispersion cube
    if edisp:
        edispcube = ctools.ctedispcube(obs)
        edispcube['incube']    = 'NONE'
        edispcube['usepnt']    = usepnt
        edispcube['ebinalg']   = 'LOG'
        edispcube['binsz']     = psf_binsz
        edispcube['nxpix']     = psf_nxpix
        edispcube['nypix']     = psf_nypix
        edispcube['enumbins']  = _enumbins
        edispcube['emin']      = emin
        edispcube['emax']      = emax
        edispcube['coordsys']  = coordsys
        edispcube['proj']      = proj
        edispcube['addbounds'] = addbounds
        edispcube['debug']     = debug
        edispcube['chatter']   = chatter
        if not usepnt:
            edispcube['xref'] = xref
            edispcube['yref'] = yref
        if log:
            edispcube.logFileOpen()
        edispcube.run()

    # Build response dictionary
    response = {}
    response['expcube'] = expcube.expcube().copy()
    response['psfcube'] = psfcube.psfcube().copy()
    response['bkgcube'] = bkgcube.bkgcube().copy()
    response['models']  = bkgcube.models().copy()
    if edisp:
        response['edispcube'] = edispcube.expcube().copy()

    # Return response cubes
    return response
