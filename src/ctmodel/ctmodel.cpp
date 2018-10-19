/***************************************************************************
 *                  ctmodel - Model cube generation tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2018 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file ctmodel.cpp
 * @brief Model cube generation tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctmodel.hpp"
#include "GTools.hpp"

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#endif

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                          "ctmodel::get_parameters()"
#define G_GET_OBS                                        "ctmodel::get_obs()"
#define G_FILL_CUBE                    "ctmodel::fill_cube(GCTAObservation*)"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Constants __________________________________________________________ */
const GEnergy g_energy_margin(1.0e-12, "TeV");


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
ctmodel::ctmodel(void) : ctobservation(CTMODEL_NAME, VERSION)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observations constructor
 *
 * param[in] obs Observation container.
 *
 * This method creates an instance of the class by copying an existing
 * observations container.
 ***************************************************************************/
ctmodel::ctmodel(const GObservations& obs) :
         ctobservation(CTMODEL_NAME, VERSION, obs)
{
    // Initialise members
    init_members();

    // Return
    return;
}



/***********************************************************************//**
 * @brief Command line constructor
 *
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 ***************************************************************************/
ctmodel::ctmodel(int argc, char *argv[]) :
         ctobservation(CTMODEL_NAME, VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Application.
 ***************************************************************************/
ctmodel::ctmodel(const ctmodel& app) : ctobservation(app)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(app);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
ctmodel::~ctmodel(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] app Application.
 * @return Application.
 ***************************************************************************/
ctmodel& ctmodel::operator=(const ctmodel& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->ctobservation::operator=(app);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(app);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear ctmodel tool
 *
 * Clears ctmodel tool.
 ***************************************************************************/
void ctmodel::clear(void)
{
    // Free members
    free_members();
    this->ctobservation::free_members();
    this->ctool::free_members();

    // Clear base class (needed to conserve tool name and version)
    this->GApplication::clear();

    // Initialise members
    this->ctool::init_members();
    this->ctobservation::init_members();
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generate the model map(s)
 *
 * This method reads the task parameters from the parfile, sets up the
 * observation container, loops over all CTA observations in the container
 * and generates a model map for each CTA observation.
 ***************************************************************************/
void ctmodel::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Extract cube properties
    extract_cube_properties();

    // Set energy dispersion flags of all CTA observations and save old
    // values in save_edisp vector
    std::vector<bool> save_edisp = set_edisp(m_obs, m_apply_edisp);

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write input model container into logger
    log_models(NORMAL, m_obs.models(), "Input model");

    // Write header
    log_header1(TERSE, "Generate model cube");

    // Loop over all observations in the container
    #pragma omp parallel
    {
        // Create a models object for this core
        GModels models(m_obs.models());

        #pragma omp for schedule(static,1)
        for (int i = 0; i < m_obs.size(); ++i) {

            // Get CTA observation
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

            // Skip observation if it's not CTA
            if (obs == NULL) {
                log_header3(TERSE, get_obs_header(m_obs[i]));
                std::string msg = " Skipping "+m_obs[i]->instrument()+
                " observation";
                log_string(NORMAL, msg);
            }

            // Fill cube and leave loop if we are binned mode (meaning we
            // only have one binned observation)
            else if (m_binned) {
                fill_cube(obs, models);
                i = m_obs.size();
            }

            // Skip observation if we have a binned observation
            else if (obs->eventtype() == "CountsCube") {
                log_header3(TERSE, get_obs_header(m_obs[i]));
                std::string msg = " Skipping binned "+m_obs[i]->instrument()+
                " observation";
                log_string(NORMAL, msg);
            }

            // Otherwise, everything seems to be fine, so fill the cube from obs
            else {
                // Fill the cube
                fill_cube(obs, models);

                // Dispose events to free memory if event file exists on disk
                if (obs->eventfile().length() > 0 && obs->eventfile().exists()) {
                    obs->dispose_events();
                }
            }

        } // endfor: looped over observations

    } // end openmp parallel section

    // Write model cube into header
    log_header1(NORMAL, "Model cube");
    log_string(NORMAL, m_cube.print());

    // Restore energy dispersion flags of all CTA observations
    restore_edisp(m_obs, save_edisp);

    // Optionally publish model cube
    if (m_publish) {
        publish();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save model cube
 *
 * Saves the model cube into a FITS file specified using the "outfile"
 * task parameter.
 ***************************************************************************/
void ctmodel::save(void)
{
    // Write header
    log_header1(TERSE, "Save model cube");

    // Save only if filename is valid and cube has some size
    if ((*this)["outcube"].is_valid() && m_cube.size() > 0) {

        // Get model cube filename
        m_outcube = (*this)["outcube"].filename();

        // Save cube
        m_cube.save(m_outcube, clobber());
    }

    // Write into logger what has been done
    std::string fname = (m_outcube.is_empty()) ? "NONE" : m_outcube.url();
    if (m_cube.size() == 0) {
        fname.append(" (cube is empty, no file created)");
    }
    log_value(NORMAL, "Model cube file", fname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Publish model cube
 *
 * @param[in] name Model cube name.
 ***************************************************************************/
void ctmodel::publish(const std::string& name)
{
    // Write header into logger
    log_header1(TERSE, "Publish model cube");

    // Set default name is user name is empty
    std::string user_name(name);
    if (user_name.empty()) {
        user_name = CTMODEL_NAME;
    }

    // Log filename
    log_value(NORMAL, "Model cube name", user_name);

    // Publish model cube
    m_cube.counts().publish(user_name);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set model cube
 *
 * @param[in] cube Model cube.
 *
 * Set model cube and set all cube bins to zero.
 ***************************************************************************/
void ctmodel::cube(const GCTAEventCube& cube)
{
    // Set cube
    m_cube = cube;

    // Set all cube bins to zero
    for (int i = 0; i < m_cube.size(); ++i) {
        m_cube[i]->counts(0.0);
    }

    // Signal that cube has been set
    m_has_cube = true;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void ctmodel::init_members(void)
{
    // Initialise members
    m_outcube.clear();
    m_apply_edisp = false;
    m_publish     = false;
    m_chatter     = static_cast<GChatter>(2);

    // Initialise protected members
    m_cube.clear();
    m_gti.clear();
    m_has_cube    = false;
    m_append_cube = false;
    m_binned      = false;
    m_dir.clear();
    m_solidangle.clear();
    m_energy.clear();
    m_ewidth.clear();
    m_time.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctmodel::copy_members(const ctmodel& app)
{
    // Copy attributes
    m_outcube     = app.m_outcube;
    m_apply_edisp = app.m_apply_edisp;
    m_publish     = app.m_publish;
    m_chatter     = app.m_chatter;

    // Copy protected members
    m_cube        = app.m_cube;
    m_gti         = app.m_gti;
    m_has_cube    = app.m_has_cube;
    m_append_cube = app.m_append_cube;
    m_binned      = app.m_binned;
    m_dir         = app.m_dir;
    m_solidangle  = app.m_solidangle;
    m_energy      = app.m_energy;
    m_ewidth      = app.m_ewidth;
    m_time        = app.m_time;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctmodel::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. The parameters are read in the correct order.
 ***************************************************************************/
void ctmodel::get_parameters(void)
{
    // Reset cube append flag
    m_append_cube = false;

    // If there are no observations in container then load them via user
    // parameters.
    if (m_obs.size() == 0) {
        get_obs();
    }

    // ... otherwise add response information and energy boundaries in case
    // that they are missing
    else {
        setup_observations(m_obs);
    }

    // If we have now excactly one CTA observation (but no cube has yet been
    // appended to the observation) then check whether this observation
    // is a binned observation, and if yes, extract the counts cube for
    // model generation
    if ((m_obs.size() == 1) && (m_append_cube == false)) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

        // Continue only if observation is a CTA observation
        if (obs != NULL) {

            // Check for binned observation
            if (obs->eventtype() == "CountsCube") {

                // Set cube from binned observation
                GCTAEventCube* evtcube = dynamic_cast<GCTAEventCube*>
                                         (const_cast<GEvents*>(obs->events()));
                cube(*evtcube);

                // Signal that cube has been set
                m_has_cube = true;

                // Signal that we are in binned mode
                m_binned = true;

            } // endif: observation was binned

        } // endif: observation was CTA

    } // endif: had exactly one observation

    // Read model definition file if required
    if (m_obs.models().size() == 0) {

        // Get model filename
        std::string inmodel = (*this)["inmodel"].filename();

        // Load models from file
        m_obs.models(inmodel);

    } // endif: there were no models

    // Get energy dispersion flag parameters
    m_apply_edisp = (*this)["edisp"].boolean();

    // If we do not have yet a counts cube for model computation then check
    // whether we should read it from the "incube" parameter or whether we
    // should create it from scratch using the task parameters
    if (!m_has_cube) {

        // If the cube filename is valid the load the cube and set all cube
        // bins to zero
        if ((*this)["incube"].is_valid()) {

            // Get cube definition filename
            std::string incube = (*this)["incube"].filename();

            // Load cube from given file
            m_cube.load(incube);

            // Set all cube bins to zero
            for (int i = 0; i < m_cube.size(); ++i) {
                m_cube[i]->counts(0.0);
            }

        } // endif: cube filename was valid

        // ... otherwise create a cube from the user parameters
        else {
            m_cube = create_cube(m_obs);
        }

        // Signal that cube has been set
        m_has_cube = true;

    } // endif: we had no cube yet

    // Get remaining parameters
    m_publish = (*this)["publish"].boolean();
    m_chatter = static_cast<GChatter>((*this)["chatter"].integer());

    // If needed later, query output filename now
    if (read_ahead()) {
        (*this)["outcube"].query();
    }

    // If cube should be appended to first observation then do that now.
    // This is a kluge that makes sure that the cube is passed as part
    // of the observation in case that a cube response is used. The kluge
    // is needed because the GCTACubeSourceDiffuse::set method needs to
    // get the full event cube from the observation. It is also at this
    // step that the GTI, which may just be a dummy GTI when create_cube()
    // has been used, will be set.
    if (m_append_cube) {

        //TODO: Check that energy boundaries are compatible

        // Attach GTI of observations to model cube
        m_cube.gti(m_obs[0]->events()->gti());

        // Attach model cube to observations
        m_obs[0]->events(m_cube);

    } // endif: cube was scheduled for appending

    // Set number of OpenMP threads
    #ifdef _OPENMP
    int nthreads = (*this)["nthreads"].integer();
    if (nthreads > 0) {
        omp_set_num_threads(nthreads);
    }
    #endif

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get observation container
 *
 * Get an observation container according to the user parameters. The method
 * supports loading of a individual FITS file or an observation definition
 * file in XML format.
 *
 * If the input filename is empty, the method checks for the existence of the
 * "expcube", "psfcube" and "bkgcube" parameters. If file names have been
 * specified, the method loads the files and creates a dummy events cube that
 * is appended to the observation container.
 *
 * If no file names are specified for the "expcube", "psfcube" or "bkgcube"
 * parameters, the method reads the necessary parameters to build a CTA
 * observation from scratch.
 *
 * The method sets m_append_cube = true and m_binned = true in case that
 * a stacked observation is requested (as detected by the presence of the
 * "expcube", "psfcube", and "bkgcube" parameters). In that case, it appended
 * a dummy event cube to the observation.
 *
 * @todo Support stacked energy dispersion
 ***************************************************************************/
void ctmodel::get_obs(void)
{
    // If no observation definition file has been specified then read all
    // parameters that are necessary to create an observation from scratch
    if ((*this)["inobs"].is_undefined()) {

        // If the filenames are valid then build an observation from cube
        // response information
        if ((*this)["expcube"].is_valid() &&
            (*this)["psfcube"].is_valid() &&
            (*this)["bkgcube"].is_valid()) {

            // Get response cube filenames
            GFilename edispcube;
            bool      query_edisp = (*this)["edisp"].boolean();
            GFilename expcube     = (*this)["expcube"].filename();
            GFilename psfcube     = (*this)["psfcube"].filename();
                if (query_edisp && (*this)["edispcube"].is_valid()) {
                edispcube = (*this)["edispcube"].filename();
            }
            GFilename bkgcube     = (*this)["bkgcube"].filename();

            // Get exposure, PSF and background cubes
            GCTACubeExposure   exposure(expcube);
            GCTACubePsf        psf(psfcube);
            GCTACubeBackground background(bkgcube);

            // Create energy boundaries
            GEbounds ebounds = create_ebounds();

            // Create dummy sky map cube
            GSkyMap map("CAR","GAL",0.0,0.0,1.0,1.0,1,1,ebounds.size());

            // Create event cube
            GCTAEventCube cube(map, ebounds, exposure.gti());

            // Create CTA observation
            GCTAObservation cta;
            cta.events(cube);

            // If querying of energy dispersion cube is requested then
            // query it now
            if (query_edisp) {

                // If filename is valid then use energy dispersion ...
                if ((*this)["edispcube"].is_valid()) {

                    // Load energy dispersion cube
                    GCTACubeEdisp edisp(edispcube);

                    // Set response with all four cubes
                    cta.response(exposure, psf, edisp, background);

                } // endif: energy dispersion cube was provided

                // ... otherwise throw an exception
                else {
                    std::string msg = "Energy dispersion requested but no "
                                      "energy dispersion cube was specified.";
                    throw GException::invalid_value(G_GET_OBS, msg);
                }

            } // endif: energy dispersion needed

            // ... otherwise work without energy dispersion
            else {
                cta.response(exposure, psf, background);
            }

            // Append observation to container
            m_obs.append(cta);

            // Signal that we are in binned mode
            m_binned = true;

            // Signal that we appended a cube
            m_append_cube = true;

        } // endif: cube response information was available

        // ... otherwise build an observation from IRF response information
        else {

            // Create CTA observation
            GCTAObservation cta = create_cta_obs();

            // Set response
            set_obs_response(&cta);

            // Append observation to container
            m_obs.append(cta);

        }

    } // endif: filename was not valid

    // ... otherwise we have a file name
    else {

        // Get the filename from the input parameters
        std::string filename = (*this)["inobs"].filename();

        // If file is a FITS file then create an empty CTA observation
        // and load file into observation
        if (GFilename(filename).is_fits()) {

            // Allocate empty CTA observation
            GCTAObservation cta;

            // Load data
            cta.load(filename);

            // Set response
            set_obs_response(&cta);

            // Append observation to container
            m_obs.append(cta);

            // Signal that no XML file should be used for storage
            m_use_xml = false;

        }

        // ... otherwise load file into observation container
        else {

            // Load observations from XML file
            m_obs.load(filename);

            // For all observations that have no response, set the response
            // from the task parameters
            set_response(m_obs);

            // Set observation boundary parameters (emin, emax, rad)
            set_obs_bounds(m_obs);

            // Signal that XML file should be used for storage
            m_use_xml = true;

        } // endelse: file was an XML file

    }

    // Return
    return;

}


/***********************************************************************//**
 * @brief Extract cube properties in data members
 ***************************************************************************/
void ctmodel::extract_cube_properties(void)
{
    // Clear cube properties
    m_dir.clear();
    m_solidangle.clear();
    m_energy.clear();
    m_ewidth.clear();

    // Get number of spatial pixels and energy layers in counts cube
    int npix   = m_cube.npix();
    int nebins = m_cube.ebins();

    // Reserve capacity in vectors
    m_dir.reserve(npix);
    m_solidangle.reserve(npix);
    m_energy.reserve(nebins);
    m_ewidth.reserve(nebins);

    // Loop over all the spatial bins and extract spatial cube properties
    for (int i = 0; i < npix; ++i) {

        // Get cube bin
        GCTAEventBin* bin = m_cube[i];

        // Extract sky direction and solid angle
        m_dir.push_back(bin->dir());
        m_solidangle.push_back(bin->solidangle());

    } // endfor: looped over spatial bins

    // Loop over all the spectral bins and extract spectral cube properties
    for (int iebin = 0, ibin = 0; iebin < nebins; ++iebin, ibin += npix) {

        // Get cube bin
        GCTAEventBin* bin = m_cube[ibin];

        // Extract energy and energy width
        m_energy.push_back(bin->energy());
        m_ewidth.push_back(bin->ewidth());

    } // endfor: looped over spectral bins

    // Set time
    m_time = m_cube[0]->time();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill model into model cube
 *
 * @param[in] obs CTA observation.
 *
 * Adds the expected number of events for a given observation to the events
 * that are already found in the model cube. The method also updates the
 * GTI of the model cube so that cube GTI is a list of the GTIs of all
 * observations that were used to generate the model cube.
 ***************************************************************************/
void ctmodel::fill_cube(const GCTAObservation* obs, GModels& models)
{
    // Get number of spatial and spectral bins in counts cube
    int npix   = m_cube.npix();
    int nebins = m_cube.ebins();

    // Get reference to the pointing for DETX and DETY computation, if
    // required
    const GCTAPointing& pnt = obs->pointing();

    // Get references to GTI and energy boundaries for the event list
    const GGti&     gti         = obs->events()->gti();
    const GEbounds& obs_ebounds = obs->ebounds();

    // Get cube energy boundaries
    const GEbounds& cube_ebounds = m_cube.ebounds();

    // Get counts cube usage flags
    std::vector<bool> usage = cube_layer_usage(cube_ebounds, obs_ebounds);

    // Initialise empty, invalid RoI
    GCTARoi roi;

    // Retrieve RoI in case we have an unbinned observation
    if (obs->eventtype() == "EventList") {
        roi = obs->roi();
    }

    // Initialise statistics
    double sum              = 0.0;
    int    num_outside_ebds = 0;
    int    num_outside_roi  = 0;

    // Extract copy of models from observation container
    GModels ex_models = trim_models(models, roi);

    // Get pointer to event cube pixels
    double* pixels = const_cast<double*>(m_cube.counts().pixels());

    // Initialise event bin
    GCTAEventBin bin;

    // Set event bin attributes that are constant for a counts cube
    bin.time(m_time);
    bin.ontime(obs->ontime());
    bin.weight(1.0);

    // Loop over all the spatial bins
    for (int i = 0; i < npix; ++i) {

        // If RoI is valid then skip bin if it is outside the RoI of the
        // observation
        if (roi.is_valid() && !roi.contains(m_dir[i])) {
            num_outside_roi += nebins;
            continue;
        }

        // Get instrument direction for spatial pixel and recompute DETX and
        // DETY based on the pointing of the observation.
        GCTAInstDir dir = pnt.instdir(m_dir[i].dir());

        // Set instrument direction, solid angle and pixel index of bin
        bin.dir(dir);
        bin.solidangle(m_solidangle[i]);
        bin.ipix(i);

        // Loop over all of the energy bins of the cube
        for (int iebin = 0, ibin = i; iebin < nebins; ++iebin, ibin += npix) {

            // Skip bin if the corresponding counts cube energy bin is not fully
            // contained in the event list energy range. This avoids having
            // partially filled bins.
            if (!usage[iebin]) {
                num_outside_ebds++;
                continue;
            }

            // Set energy, energy width, energy index and weight of bin
            bin.energy(m_energy[iebin]);
            bin.ewidth(m_ewidth[iebin]);
            bin.ieng(iebin);
            //bin.weight(m_cube[ibin]->weight()); // Bad residuals for H.E.S.S. Crab

            // Compute model value for cube bin
            double model = models.eval(bin, *obs) * bin.size();

            // Sum model
            sum += model;

            // Store value
            #pragma omp critical(ctmodel_fill_cube)
            pixels[ibin] += model;

        } // endfor: looped over all spatial bins

    } // endfor: looped over all cube bins

    // Re-append the models to the main model
    models.extend(ex_models);

    // Update GTIs
    #pragma omp critical(ctmodel_fill_cube)
    {
        // Append GTIs of observation to list of GTIs
        m_gti.extend(gti);

        // Update GTIs
        m_cube.gti(m_gti);
    }

    // Log filling results
    #pragma omp critical(ctmodel_fill_cube)
    {
        log_header3(TERSE, get_obs_header(obs));
        log_value(NORMAL, "Model events in cube", sum);
        log_value(NORMAL, "Bins outside energy range", num_outside_ebds);
        log_value(NORMAL, "Bins outside RoI", num_outside_roi);
    }

    // Write model cube into header
    if (m_chatter >= EXPLICIT) {
        log_header2(EXPLICIT, "Model cube");
        log_string(EXPLICIT, m_cube.print(m_chatter));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Find the models falling inside a defined region of interest.
 *
 * @param[in,out] all_models Model container to be trimmed.
 * @param[in] roi Observation region of interest.
 * @return New model container containing only models inside RoI.
 *
 * Note that a buffer is added to the observation region to ensure that
 * point sources (which have radius=0) are appropriately included if they
 * fall near the edge of the observation.
 ***************************************************************************/
GModels ctmodel::trim_models(GModels& all_models, const GCTARoi& roi)
{
    // Create the model container containing only models in RoI
    GModels excluded_models;

    // Do nothing if roi is not valid. This will be the case for binned data
    // sets as they have the default roi radius of 0. In this case we would
    // remove ALL models, which is obviously incorrect. So we keep all models
    // when filling based on a counts cube.
    if (!roi.is_valid()) {
        return excluded_models;
    }

    // Remove all models that dont overlap with the region of interest. Note
    // that an extra factor is used since point sources have regions of
    // radius 0, which is a problem when they fall just barely outside the
    // RoI. This ensures point sources are appropriately included.
    GSkyRegionCircle obsreg(roi.centre().dir(), roi.radius()+0.5);

    // Loop over the models in the passed container
    for (int i = 0; i < all_models.size(); ++i) {

        // Cast the model to a GModelSky object
        GModelSky* model = dynamic_cast<GModelSky*>(all_models[i]);

        // If model is not a GModelSky model than skip
        if (model == NULL) {
            continue;
        }

        // Otherwise, exclude the model if it isn't a diffuse cube and it
        // doesn't overlap with the observation.
        if (!(model->spatial()->code() == GMODEL_SPATIAL_DIFFUSE) &&
            (!model->spatial()->region()->overlaps(obsreg))) {

            // Append model to the excluded models
            excluded_models.append(*model);

            // Remove model
            all_models.remove(i);

            // Decrement to prevent skipping models
            i--;
        }

    } // endfor: looped over all models

    // Return excluded models
    return excluded_models;
}
