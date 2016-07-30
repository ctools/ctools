# ==========================================================================
# CTA observation handling support functions
#
# Copyright (C) 2011-2016 Juergen Knoedlseder
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


# ===================== #
# Simulate observations #
# ===================== #
def sim(obs, log=False, debug=False, chatter=2, edisp=False, seed=0,
        emin=None, emax=None, nbins=0,
        binsz=0.05, npix=200, proj='TAN', coord='GAL'):
    """
    Simulate events for all observations in the container

    Simulate events for all observations using ctobssim. If the number of
    energy bins is positive, the events are then binned in a counts cube
    using ctbin. If multiple observations are simulated, the counts cube is
    a stacked cube and the corresponding response cubes are computed using
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

        # If energy boundaries are not given then determine common energy
        # boundaries for all observations
        if emin == None or emax == None:
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

            # First get the cube. We need this because ctbin builds a cube on
            # the fly. We should avoid this and make the cube part of ctbin
            cube = bin.cube()

            # Get xref, yref, coordsys and proj from counts cube
            _xref     = cube.counts().projection().crval(0)
            _yref     = cube.counts().projection().crval(1)
            _coordsys = cube.counts().projection().coordsys()
            _proj     = cube.counts().projection().code()

            # Kluge: map "EQU" to "CEL" coordinate system (not sure why the
            # GSkyProjection::coordsys returns "EQU" instead of "CEL"
            if _coordsys == 'EQU':
                _coordsys = 'CEL'

            # Get stacked response
            response = get_stacked_response(sim.obs(), _xref, _yref,
                                            binsz=binsz, nxpix=npix, nypix=npix,
                                            emin=emin, emax=emax, enumbins=nbins,
                                            edisp=edisp,
                                            coordsys=_coordsys, proj=_proj,
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


# ================ #
# Fit observations #
# ================ #
def fit(obs, log=False, debug=False, chatter=2, edisp=False):
    """
    Perform maximum likelihood fitting of observations in the container

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    log : bool, optional
        Create log file(s)
    debug : bool, optional
        Create console dump?
    chatter : int, optional
        Chatter level
    edisp : bool, optional
        Apply energy dispersion?

    Returns
    -------
    like : `~ctools.ctlike`
        ctlike application
    """
    # Allocate ctlike application and set parameters
    like = ctools.ctlike(obs)
    like['debug']   = debug
    like['chatter'] = chatter
    like['edisp']   = edisp

    # Optionally open the log file
    if log:
        like.logFileOpen()

    # Run ctlike application.
    like.run()

    # Return observations
    return like


# ============================================================== #
# Fit observations and determine errors using likelihood profile #
# ============================================================== #
def cterror(obs, srcname, log=False, debug=False, chatter=2):
    """
    Fit observations and determine errors using likelihood profile

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    srcname : str
        Source name
    log : bool, optional
        Create log file(s)
    debug : bool, optional
        Create console dump?
    chatter : int, optional
        Chatter level

    Returns
    -------
    error : `~ctools.cterror`
        cterror application
    """
    # Allocate cterror application and set parameters
    error = ctools.cterror(obs)
    error['srcname'] = srcname
    error['debug']   = debug
    error['chatter'] = chatter

    # Optionally open the log file
    if log:
        error.logFileOpen()

    # Run cterror application.
    error.run()

    # Return observations
    return error


# ================= #
# Create counts map #
# ================= #
def cntmap(obs, proj='TAN', coord='GAL', xval=0.0, yval=0.0, \
           binsz=0.05, nxpix=200, nypix=200, \
           outname='cntmap.fits'):
    """
    Creates a counts map by combining the events of all observations

    The counts map will be a summed map over all energies

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container without events
    proj : str, optional
        Projection for binned simulation
    coord : str, optional
        Coordinate system for binned simulation
    xval : float, optional
        Counts map centre of X direction (deg)
    yval : float, optional
        Counts map centre of Y direction (deg)
    binsz : float, optional
        Pixel size for binned simulation (deg/pixel)
    nxpix : int, optional
        Number of pixels in X direction
    nypix : int, optional
        Number of pixels in Y direction
    outname : str, optional
        Counts map file name

    Returns
    -------
    map : `~gammalib.GSkyMap`
        Counts map
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
def modmap(obs, eref=0.1, proj='TAN', coord='GAL', xval=0.0, yval=0.0, \
           binsz=0.05, nxpix=200, nypix=200, \
           outname='modmap.fits'):
    """
    Make model map for a given reference energy by combining all observations

    The model map will be evaluated for a given reference energy 'eref' and will
    be given in units of [counts/(sr MeV s)].

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container without events
    eref : float, optional
        Reference energy for which model is created (TeV)
    proj : str, optional
        Projection for binned simulation
    coord : str, optional
        Coordinate system for binned simulation
    xval : float, optional
        Counts map centre of X direction (deg)
    yval : float, optional
        Counts map centre of Y direction (deg)
    binsz : float, optional
        Pixel size for binned simulation (deg/pixel)
    nxpix : int, optional
        Number of pixels in X direction
    nypix : int, optional
        Number of pixels in Y direction
    outname : str, optional
        Model map file name

    Returns
    -------
    map : `~gammalib.GSkyMap`
        Model map
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


# ======================= #
# Set one CTA observation #
# ======================= #
def set_obs(pntdir, tstart=0.0, duration=1800.0, deadc=0.95, \
            emin=0.1, emax=100.0, rad=5.0, \
            irf='South_50h', caldb='prod2', id='000000'):
    """
    Set a single CTA observation
    
    The function sets a single CTA observation containing an empty CTA
    event list. By looping over this function you can add CTA observations
    to the observation container.

    Parameters
    ----------
    pntdir : `~gammalib.GSkyDir`
        Pointing direction
    tstart : float, optional
        Start time (seconds)
    duration : float, optional
        Duration of observation (seconds)
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
    id : str, optional
        Observation identifier

    Returns
    -------
    obs : `~gammalib.GCTAObservation`
        CTA observation
    """
    # Allocate CTA observation
    obs_cta = gammalib.GCTAObservation()

    # Set calibration database
    db = gammalib.GCaldb()
    if (gammalib.dir_exists(caldb)):
        db.rootdir(caldb)
    else:
        db.open('cta', caldb)

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
    ebounds = gammalib.GEbounds(gammalib.GEnergy(emin, 'TeV'),
                                gammalib.GEnergy(emax, 'TeV'))

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
        Start time (seconds)
    duration : float, optional
        Duration of observation (seconds)
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
        id = '%6.6d' % obs_id

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


# ======================== #
# Set model for TS fitting #
# ======================== #
def set_ts_model(models, srcname, ra=None, dec=None, fitspat=False, fitspec=False):
    """
    Set model for TS fitting

    Parameters
    ----------
    models : `~gammalib.GModels`
        Input model container
    srcname : str
        Test source name
    ra : float, optional
        Right Ascension of test source (deg)
    dec : float, optional
        Declination of test source (deg)
    fitspat : bool, optional
        Fit spatial parameter?
    fitspec : bool, optional
        Fit spectral parameters?

    Returns
    -------
    model : `~gammalib.GModels`
        Model container for TS fitting
    """
    # Create a clone of the input model
    outmodels = models.copy()

    # Disable TS computation for all model components (will enable the
    # test source later)
    for model in outmodels:
        model.tscalc(False)

    # Get source model and enable TS computation
    model = outmodels[srcname]
    model.tscalc(True)

    # If source model has no "Prefactor" parameter then raise an exception
    if not model.has_par('Prefactor'):
        msg = ('Model "%s" has no parameter "Prefactor". Only spectral '
               'models with a "Prefactor" parameter are supported.' %
               srcname)
        raise RuntimeError(msg)

    # Set position of test source
    if ra != None and dec != None:
        if model.has_par('RA') and model.has_par('DEC'):
            model['RA'].value(ra)
            model['DEC'].value(dec)

    # Set possible spatial and spectral parameters
    spatial  = ['RA', 'DEC', 'Sigma', 'Radius', 'Width', 'PA',
                'MinorRadius', 'MajorRadius']
    spectral = ['Index', 'Index1', 'Index2', 'BreakEnergy', 'CutoffEnergy',
                'InverseCutoffEnergy']

    # Fit or fix spatial parameters
    for par in spatial:
        if model.has_par(par):
            if fitspat:
                model[par].free()
            else:
                model[par].fix()

    # Fit or fix spectral parameters
    for par in spectral:
        if model.has_par(par):
            if fitspec:
                model[par].free()
            else:
                model[par].fix()

    # Return model container
    return outmodels


# ==================== #
# Get stacked response #
# ==================== #
def get_stacked_response(obs, xref, yref, binsz=0.05, nxpix=200, nypix=200,
                         emin=0.1, emax=100.0, enumbins=20, edisp=False,
                         coordsys='GAL', proj='TAN',
                         log=False, debug=False, chatter=2):
    """
    Get stacked response cubes

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
    # Compute spatial binning for point spread function and
    # energy dispersion cubes
    psf_binsz = 10.0 * binsz
    psf_nxpix = min(nxpix // 10, 2)  # Make sure result is int
    psf_nypix = min(nypix // 10, 2)  # Make sure result is int

    # Create exposure cube
    expcube = ctools.ctexpcube(obs)
    expcube['incube']   = 'NONE'
    expcube['usepnt']   = False
    expcube['ebinalg']  = 'LOG'
    expcube['xref']     = xref
    expcube['yref']     = yref
    expcube['binsz']    = binsz
    expcube['nxpix']    = nxpix
    expcube['nypix']    = nypix
    expcube['enumbins'] = enumbins
    expcube['emin']     = emin
    expcube['emax']     = emax
    expcube['coordsys'] = coordsys
    expcube['proj']     = proj
    if log:
        expcube.logFileOpen()
    expcube.run()

    # Create point spread function cube
    psfcube = ctools.ctpsfcube(obs)
    psfcube['incube']   = 'NONE'
    psfcube['usepnt']   = False
    psfcube['ebinalg']  = 'LOG'
    psfcube['xref']     = xref
    psfcube['yref']     = yref
    psfcube['binsz']    = psf_binsz
    psfcube['nxpix']    = psf_nxpix
    psfcube['nypix']    = psf_nypix
    psfcube['enumbins'] = enumbins
    psfcube['emin']     = emin
    psfcube['emax']     = emax
    psfcube['coordsys'] = coordsys
    psfcube['proj']     = proj
    if log:
        psfcube.logFileOpen()
    psfcube.run()

    # Create background cube
    bkgcube = ctools.ctbkgcube(obs)
    bkgcube['incube']   = 'NONE'
    bkgcube['usepnt']   = False
    bkgcube['ebinalg']  = 'LOG'
    bkgcube['xref']     = xref
    bkgcube['yref']     = yref
    bkgcube['binsz']    = binsz
    bkgcube['nxpix']    = nxpix
    bkgcube['nypix']    = nypix
    bkgcube['enumbins'] = enumbins
    bkgcube['emin']     = emin
    bkgcube['emax']     = emax
    bkgcube['coordsys'] = coordsys
    bkgcube['proj']     = proj
    if log:
        bkgcube.logFileOpen()
    bkgcube.run()

    # If energy dispersion is requested then create energy dispersion cube
    if edisp:
        edispcube = ctools.ctedispcube(obs)
        edispcube['incube']   = 'NONE'
        edispcube['usepnt']   = False
        edispcube['ebinalg']  = 'LOG'
        edispcube['xref']     = xref
        edispcube['yref']     = yref
        edispcube['binsz']    = psf_binsz
        edispcube['nxpix']    = psf_nxpix
        edispcube['nypix']    = psf_nypix
        edispcube['enumbins'] = enumbins
        edispcube['emin']     = emin
        edispcube['emax']     = emax
        edispcube['coordsys'] = coordsys
        edispcube['proj']     = proj
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
