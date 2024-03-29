# ==========================================================================
# Utility functions for observation handling
#
# Copyright (C) 2011-2020 Juergen Knoedlseder
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
import cscripts


# ===================== #
# Simulate observations #
# ===================== #
def sim(obs, log=False, debug=False, chatter=2, edisp=False, seed=0,
        emin=None, emax=None, nbins=0, onsrc=None, onrad=0.2, addbounds=False,
        binsz=0.05, npix=200, proj='TAN', coord='GAL', nthreads=0):
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
    onsrc : str, optional
        Name of source for On region (None if no On/Off obs. is used)
    onrad : float, optional
        Radius for On region (deg)
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
    nthreads : str, optional
        Number of parallel processes for On/Off spectra computation (0=all available CPUs)

    Returns
    -------
    obs : `~gammalib.GObservations`
        Observation container filled with simulated events
    """
    # Allocate ctobssim application and set parameters
    obssim = ctools.ctobssim(obs)
    obssim['seed']     = seed
    obssim['edisp']    = edisp
    obssim['nthreads'] = nthreads
    obssim['chatter']  = chatter
    obssim['debug']    = debug

    # Optionally open the log file
    if log:
        obssim.logFileOpen()

    # Run ctobssim application. This will loop over all observations in the
    # container and simulation the events for each observation. Note that
    # events are not added together, they still apply to each observation
    # separately.
    obssim.run()

    # Binned option?
    if nbins > 0:

        # If energy boundaries are not given then determine the minimum and
        # the maximum energies from all observations and use these values
        # as energy boundaries. The energy boundaries are given in TeV.
        if emin == None or emax == None:
            emin = 1.0e30
            emax = 0.0
            for run in obssim.obs():
                emin = min(run.events().ebounds().emin().TeV(), emin)
                emax = max(run.events().ebounds().emax().TeV(), emax)

        # If a On source is specified then create On/Off observations
        if onsrc != None:

            # Extract source position from model
            model = obssim.obs().models()[onsrc]
            ra    = model['RA'].value()
            dec   = model['DEC'].value()

            # Allocate csphagen application and set parameters
            phagen = cscripts.csphagen(obssim.obs())
            phagen['inmodel']     = 'NONE'
            phagen['srcname']     = onsrc
            phagen['ebinalg']     = 'LOG'
            phagen['emin']        = emin
            phagen['emax']        = emax
            phagen['enumbins']    = nbins
            phagen['srcshape']    = 'CIRCLE'
            phagen['coordsys']    = 'CEL'
            phagen['ra']          = ra
            phagen['dec']         = dec
            phagen['rad']         = onrad
            phagen['stack']       = False
            phagen['inexclusion'] = 'NONE'
            phagen['bkgmethod']   = 'REFLECTED'
            phagen['nthreads']    = nthreads

            # Optionally open the log file
            if log:
                phagen.logFileOpen()

            # Run csphagen application
            phagen.run()

            # Make a deep copy of the observation that will be returned
            # (the csphagen object will go out of scope one the function is
            # left)
            obs = phagen.obs().copy()

        # ... otherwise use binned observations
        else:

            # Allocate ctbin application and set parameters
            binning = ctools.ctbin(obssim.obs())
            binning['ebinalg']  = 'LOG'
            binning['emin']     = emin
            binning['emax']     = emax
            binning['enumbins'] = nbins
            binning['usepnt']   = True # Use pointing for map centre
            binning['nxpix']    = npix
            binning['nypix']    = npix
            binning['binsz']    = binsz
            binning['coordsys'] = coord
            binning['proj']     = proj
            binning['nthreads'] = nthreads
            binning['chatter']  = chatter
            binning['debug']    = debug

            # Optionally open the log file
            if log:
                binning.logFileOpen()

            # Run ctbin application. This will loop over all observations in
            # the container and bin the events in counts maps
            binning.run()

            # If we have multiple input observations then create stacked response
            # cubes and append them to the observation
            if len(obssim.obs()) > 1:

                # Get counts cube. The counts cube is needed to obtained a
                # proper background cube.
                cntcube = binning.cube()

                # Get stacked response
                response = get_stacked_response(obssim.obs(), cntcube,
                                                edisp=edisp,
                                                addbounds=addbounds,
                                                log=log, debug=debug,
                                                chatter=chatter)

                # Set stacked response
                if edisp:
                    binning.obs()[0].response(response['expcube'],
                                              response['psfcube'],
                                              response['edispcube'],
                                              response['bkgcube'])
                else:
                    binning.obs()[0].response(response['expcube'],
                                              response['psfcube'],
                                              response['bkgcube'])

                # Set new models
                binning.obs().models(response['models'])

            # Make a deep copy of the observation that will be returned
            # (the ctbin object will go out of scope one the function is
            # left)
            obs = binning.obs().copy()

    else:

        # Make a deep copy of the observation that will be returned
        # (the ctobssim object will go out of scope one the function is
        # left)
        obs = obssim.obs().copy()

    # Delete the simulation
    del obssim

    # Return observation container
    return obs


# ======================= #
# Set one CTA observation #
# ======================= #
def set_obs(pntdir, tstart=0.0, duration=1800.0, deadc=0.98, \
            emin=0.1, emax=100.0, rad=5.0, \
            instrument='CTA', irf='South_50h', caldb='prod2', \
            obsid='000000', mjdref=51544.5):
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
    instrument : str, optional
        Name of Cherenkov Telescope
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

    # Set mission
    mission = gammalib.toupper(instrument)
    obs.instrument(mission)

    # Set CTA calibration database
    db = gammalib.GCaldb()
    if (gammalib.dir_exists(caldb)):
        db.rootdir(caldb)
    else:
        db.open(mission, caldb)

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
    ref = gammalib.GTimeReference(mjdref,'s','TT','LOCAL')
    gti = gammalib.GGti(ref)
    gti.append(gammalib.GTime(tstart, ref),
               gammalib.GTime(tstart+duration, ref))

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

    # ... otherwise raise an exception since we have an unknown pattern
    else:
        msg = 'Observation pattern "%s" not recognized.' % (pattern)
        raise RuntimeError(msg)

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
def get_stacked_response(obs, cntcube, edisp=False, addbounds=False,
                         log=False, debug=False, chatter=2):
    """
    Get stacked response cubes

    The number of energies bins are set to at least 30 bins per decade, unless
    the counts cube has more energy bins per decade.

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    edisp : bool, optional
        Apply energy dispersion?
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
    # Derive binning parameters from counts cube
    xref     = cntcube.counts().projection().crval(0)
    yref     = cntcube.counts().projection().crval(1)
    binsz    = cntcube.counts().projection().cdelt(1)
    coordsys = cntcube.counts().projection().coordsys()
    proj     = cntcube.counts().projection().code()
    nxpix    = cntcube.nx()
    nypix    = cntcube.ny()
    emin     = cntcube.emin().TeV()
    emax     = cntcube.emax().TeV()
    ebins    = cntcube.ebins()

    # Translate coordinate system
    if coordsys == 'EQU':
        coordsys = 'CEL'

    # Set energy limits, with larger etrue limits for the case that energy
    # dispersion is requested
    if edisp:
        emin *= 0.5
        emax *= 1.5

    # Set number of energy bins to at least 30 per energy decade
    enumbins = int((math.log10(emax) - math.log10(emin)) * 30.0)
    if ebins > enumbins:
        enumbins = ebins

    # Compute spatial binning for point spread function and energy dispersion
    # cubes. The spatial binning is 10 times coarser than the spatial binning
    # of the exposure and background cubes. At least 2 spatial are required.
    psf_binsz = 10.0 * binsz
    psf_nxpix = max(nxpix // 10, 2)  # Make sure result is int
    psf_nypix = max(nypix // 10, 2)  # Make sure result is int

    # Create exposure cube
    expcube = ctools.ctexpcube(obs)
    expcube['incube']    = 'NONE'
    expcube['ebinalg']   = 'LOG'
    expcube['xref']      = xref
    expcube['yref']      = yref
    expcube['binsz']     = binsz
    expcube['nxpix']     = nxpix
    expcube['nypix']     = nypix
    expcube['enumbins']  = enumbins
    expcube['emin']      = emin
    expcube['emax']      = emax
    expcube['coordsys']  = coordsys
    expcube['proj']      = proj
    expcube['addbounds'] = addbounds
    expcube['debug']     = debug
    expcube['chatter']   = chatter
    if log:
        expcube.logFileOpen()
    expcube.run()

    # Create point spread function cube
    psfcube = ctools.ctpsfcube(obs)
    psfcube['incube']    = 'NONE'
    psfcube['ebinalg']   = 'LOG'
    psfcube['xref']      = xref
    psfcube['yref']      = yref
    psfcube['binsz']     = psf_binsz
    psfcube['nxpix']     = psf_nxpix
    psfcube['nypix']     = psf_nypix
    psfcube['enumbins']  = enumbins
    psfcube['emin']      = emin
    psfcube['emax']      = emax
    psfcube['coordsys']  = coordsys
    psfcube['proj']      = proj
    psfcube['addbounds'] = addbounds
    psfcube['debug']     = debug
    psfcube['chatter']   = chatter
    if log:
        psfcube.logFileOpen()
    psfcube.run()

    # Create background cube. Note that we use the same energy binning as
    # for the counts cube since no interpolation should be actually done
    # for the background cube.
    bkgcube = ctools.ctbkgcube(obs)
    bkgcube.cntcube(cntcube)
    bkgcube['debug']   = debug
    bkgcube['chatter'] = chatter
    if log:
        bkgcube.logFileOpen()
    bkgcube.run()

    # If energy dispersion is requested then create energy dispersion cube
    if edisp:
        edispcube = ctools.ctedispcube(obs)
        edispcube['incube']    = 'NONE'
        edispcube['ebinalg']   = 'LOG'
        edispcube['xref']      = xref
        edispcube['yref']      = yref
        edispcube['binsz']     = psf_binsz
        edispcube['nxpix']     = psf_nxpix
        edispcube['nypix']     = psf_nypix
        edispcube['enumbins']  = enumbins
        edispcube['emin']      = emin
        edispcube['emax']      = emax
        edispcube['coordsys']  = coordsys
        edispcube['proj']      = proj
        edispcube['addbounds'] = addbounds
        edispcube['debug']     = debug
        edispcube['chatter']   = chatter
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
        response['edispcube'] = edispcube.edispcube().copy()

    # Return response cubes
    return response


# ================================= #
# Get stacked observation container #
# ================================= #
def get_stacked_obs(cls, obs, nthreads=0):
    """
    Bin an observation and return an observation container with a single
    binned observation

    Parameters
    ----------
    cls : `~ctools.cscript`
        cscript class
    obs : `~gammalib.GObservations`
        Observation container
    nthreads : str, optional
        Number of parallel processes for On/Off spectra computation (0=all available CPUs)

    Returns
    -------
    obs : `~gammalib.GObservations`
        Observation container where the first observation is a binned observation
    """
    # Write header
    if cls._logExplicit():
        cls._log.header3('Binning events')

    # Bin events
    cntcube = ctools.ctbin(obs)
    cntcube['usepnt']   = False
    cntcube['ebinalg']  = 'LOG'
    cntcube['xref']     = cls['xref'].real()
    cntcube['yref']     = cls['yref'].real()
    cntcube['binsz']    = cls['binsz'].real()
    cntcube['nxpix']    = cls['nxpix'].integer()
    cntcube['nypix']    = cls['nypix'].integer()
    cntcube['enumbins'] = cls['enumbins'].integer()
    cntcube['emin']     = cls['emin'].real()
    cntcube['emax']     = cls['emax'].real()
    cntcube['coordsys'] = cls['coordsys'].string()
    cntcube['proj']     = cls['proj'].string()
    cntcube['nthreads'] = nthreads
    cntcube.run()

    # Store counts cube so that we can use it to build the background
    # cube
    cube = cntcube.cube()

    # Write header
    if cls._logExplicit():
        cls._log.header3('Creating stacked response')

    # Get stacked response
    response = get_stacked_response(obs, cube, edisp=cls['edisp'].boolean())

    # Retrieve a new oberservation container
    new_obs = cntcube.obs().copy()

    # Set stacked response
    if cls['edisp'].boolean():
        new_obs[0].response(response['expcube'], response['psfcube'],
                            response['edispcube'], response['bkgcube'])
    else:
        new_obs[0].response(response['expcube'], response['psfcube'],
                            response['bkgcube'])

    # Get new models
    models = response['models']

    # Set models for new oberservation container
    new_obs.models(models)

    # Return new oberservation container
    return new_obs


# ================================ #
# Get On/Off observation container #
# ================================ #
def get_onoff_obs(cls, obs, nthreads=0, ra = None, dec = None, srcname = ''):
    """
    Create On/Off observations container from given observations

    Parameters
    ----------
    cls : `~ctools.cscript`
        cscript class
    obs : `~gammalib.GObservations`
        Observation container
    nthreads : str, optional
        Number of parallel processes for On/Off spectra computation (0=all available CPUs)
    ra : float, optional
        R.A. of source region centre
    dec : float, optional
        Dec. of source region centre
    srcname : str, optional
        Source name

    Returns
    -------
    onoff_obs : `~gammalib.GObservations`
        Observation container with On/Off observations
    """
    # Write header
    if cls._logExplicit():
        cls._log.header3('Creating On/Off observations')

    # Initialise inmodel and use_model_bkg
    inmodel       = 'NONE'
    use_model_bkg = True

    # Set inmodel, srcname, and use_model_bkg if they are available
    if cls.has_par('inmodel') and obs.models().size() == 0:
        if cls['inmodel'].is_valid():
            inmodel = cls['inmodel'].value()
    if srcname == '' and cls.has_par('srcname'):
        srcname = cls['srcname'].value()
    if cls.has_par('use_model_bkg'):
        use_model_bkg = cls['use_model_bkg'].boolean()

    # Create On/Off observations
    phagen = cscripts.csphagen(obs)
    phagen['inmodel']     = inmodel
    phagen['srcname']     = srcname
    phagen['emin']        = cls['emin'].real()
    phagen['emax']        = cls['emax'].real()
    phagen['enumbins']    = cls['enumbins'].integer()
    phagen['ebinalg']     = 'LOG'
    phagen['srcshape']    = cls['srcshape'].string()
    # User has specified custom centre for the source region
    if ra != None and dec != None:
        phagen['coordsys'] = 'CEL'
        phagen['ra']       = ra
        phagen['dec']      = dec
    # Otherwise use default in class
    else:
        phagen['coordsys']    = cls['coordsys'].string()
        if cls['coordsys'].string() == 'CEL':
            phagen['ra']  = cls['xref'].real()
            phagen['dec'] = cls['yref'].real()
        elif cls['coordsys'].string() == 'GAL':
            phagen['glon'] = cls['xref'].real()
            phagen['glat'] = cls['yref'].real()
    if cls['srcshape'].string() == 'CIRCLE':
        phagen['rad'] = cls['rad'].real()
    elif cls['srcshape'].string() == 'RECT':
        phagen['width']  = cls['width'].real()
        phagen['height'] = cls['height'].real()
        phagen['posang'] = cls['posang'].real()
    phagen['bkgmethod'] = cls['bkgmethod'].string()
    if cls['bkgmethod'].string() == 'REFLECTED':
        phagen['bkgregmin'] = cls['bkgregmin'].integer()
    phagen['use_model_bkg'] = use_model_bkg
    phagen['maxoffset']     = cls['maxoffset'].real()
    phagen['stack']         = True
    phagen['etruemin']      = cls['etruemin'].real()
    phagen['etruemax']      = cls['etruemax'].real()
    phagen['etruebins']     = cls['etruebins'].integer()
    phagen['chatter']       = cls['chatter'].integer()
    phagen['clobber']       = cls['clobber'].boolean()
    phagen['debug']         = cls['debug'].boolean()
    phagen['nthreads']      = nthreads

    # Set exclusion map
    # Initialise exclusion map flag to False
    use_excl_map = False
    # Initialise inexclusion par flag to False
    use_inexclusion = False
    # Do we have valid exclusion map in memory?
    if hasattr(cls, 'exclusion_map'):
        exclusion_map = cls.exclusion_map()
        if exclusion_map is not None:
            # Set use exclusion map flag to True
            use_excl_map = True
            # Set exclusion map
            phagen.exclusion_map(exclusion_map)
    # If we do not have valid exclusion map in memory ...
    if not use_excl_map:
        # ... do we have an inxeclusion parameter?
        if cls.has_par('inexclusion'):
            # If the inexclusion parameter is valid
            if cls['inexclusion'].is_valid():
                # Set use inexclusion flag to True
                use_inexclusion = True
                # Set inexclusion parameter
                inexclusion = cls['inexclusion'].value()
                phagen['inexclusion'] = inexclusion
    # If there is no valid exclusion map in memory
    # nor valid inexclusion parameter
    if not use_excl_map and not use_inexclusion:
        # Set inexclusion for csphagen to None
        phagen['inexclusion'] = 'NONE'

    # Run csphagen
    phagen.run()

    # Clone resulting observation container
    onoff_obs = phagen.obs().copy()

    # On/Off observations are created with CSTAT as default statistic
    # We will change this to the user choice
    if cls['statistic'].string() != 'DEFAULT':
        for observation in onoff_obs:
            observation.statistic(cls['statistic'].string())

    # Return On/Off oberservation container
    return onoff_obs


# ========================================= #
# Calculate residuals from counts and model #
# ========================================= #
def residuals(cls, counts, model):
    """
    Calculate residuals given counts and models, according to algorithm
    specified by user.

    Can handle GSkyMap or GNdarray objects


    Parameters
    ----------
    cls : `~ctools.cscript`
        cscript class
    counts : `~gammalib.GSkyMap/~gammalib.GNdarray'
        Data counts
    model : `~gammalib.GSkyMap/~gammalib.GNdarray'
        Model counts

    Returns
    -------
    residuals : `~gammalib.GSkyMap/~gammalib.GNdarray'
        Residuals
    """
    # Find type of objects we are manipulating and set size to iterate later
    if counts.classname() == 'GNdarray':
        nelem = counts.size()
    elif counts.classname() == 'GSkyMap':
        nelem = counts.npix()
    else:
        msg = 'cscripts.obsutils.residuals only handles ' \
              + 'gammalib.GNdarray or gammalib.GSkyMap objects.\n'
        raise RuntimeError(msg)

    # Copy counts to initialise residuals
    residuals = counts.copy()

    # Get residual map algorithm type
    algorithm = cls['algorithm'].string()

    # Subtract
    if algorithm == 'SUB':
        residuals -= model

    # Subtract and divide by model
    elif algorithm == 'SUBDIV':
        residuals -= model
        residuals /= model

    # Subtract and divide by sqrt of model
    elif algorithm == 'SUBDIVSQRT':
        residuals -= model
        residuals /= model.sqrt()

    # Calculate significance from Li&Ma
    elif algorithm == 'SIGNIFICANCE':

        # Compute sign
        sign = (residuals - model).sign()

        # Loop over every bin
        for i in range(nelem):

            # If the model value > 0.0 do the computation as normal ...
            model_val = model[i]
            if model_val > 0.0:

                # If the data value is also > 0 then compute the
                # significance^2 and save it ...
                data_val = residuals[i]
                if data_val > 0.0:
                    log_val      = math.log(data_val / model_val)
                    residuals[i] = (data_val * log_val) + model_val - data_val
                    if residuals[i] < 0.0:   # See glitch issue #2765
                        residuals[i] = 0.0

                # ... otherwise compute the reduced value of the above
                # expression. This is necessary to avoid computing log(0).
                else:
                    residuals[i] = model_val

            # ... otherwise hard-code the significance to 0
            else:
                residuals[i] = 0.0

        # Compute significance map
        residuals *= 2.0
        residuals  = residuals.sqrt()
        residuals *= sign

    # Raise exception if algorithm is unknown
    else:
        raise TypeError('Algorithm "' + algorithm + '" not known')

    # Return
    return residuals


# ================= #
# Create counts map #
# ================= #
def create_counts_cube(cls, obs):
    """
    Create counts cube from observations

    Parameters
    ----------
    cls : `~ctools.cscript`
        cscript class
    obs : `~gammalib.GObservations`
        Observation container

    Returns
    -------
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    """
    # Initialise ctbin
    ctbin = ctools.ctbin(obs)
    
    # Set parameters
    ctbin['xref']     = cls['xref'].real()
    ctbin['yref']     = cls['yref'].real()
    ctbin['proj']     = cls['proj'].string()
    ctbin['coordsys'] = cls['coordsys'].string()
    ctbin['ebinalg']  = cls['ebinalg'].string()
    ctbin['nxpix']    = cls['nxpix'].integer()
    ctbin['nypix']    = cls['nypix'].integer()
    ctbin['binsz']    = cls['binsz'].real()
    if cls['ebinalg'].string() == 'FILE':
        ctbin['ebinfile'] = cls['ebinfile'].filename().file()
    else:
        ctbin['enumbins'] = cls['enumbins'].integer()
        ctbin['emin']     = cls['emin'].real()
        ctbin['emax']     = cls['emax'].real()
        if cls['ebinalg'].string() == 'POW':
            ctbin['ebingamma'] = cls['ebingamma'].real()
    ctbin['chatter']  = cls['chatter'].integer()
    ctbin['clobber']  = cls['clobber'].boolean()
    ctbin['debug']    = cls['debug'].boolean()

    # Run ctbin
    ctbin.run()

    # Retrieve counts cube
    cntcube = ctbin.cube().copy()

    # Return counts cube
    return cntcube
