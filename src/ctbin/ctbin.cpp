/***************************************************************************
 *                        ctbin - Event binning tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @file ctbin.cpp
 * @brief Event binning tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctbin.hpp"
#include "GTools.hpp"

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#endif

/* __ Method name definitions ____________________________________________ */
#define G_CUBE                                            "ctbin::cube(int&)"
#define G_GET_PARAMETERS                            "ctbin::get_parameters()"
#define G_FILL_CUBE                      "ctbin::fill_cube(GCTAObservation*)"
#define G_SET_WEIGHTS                  "ctbin::set_weights(GCTAObservation*)"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty event binning tool.
 ***************************************************************************/
ctbin::ctbin(void) : ctobservation(CTBIN_NAME, VERSION)
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
 * Constructs event binning tool from an observation container.
 ***************************************************************************/
ctbin::ctbin(const GObservations& obs) :
       ctobservation(CTBIN_NAME, VERSION, obs)
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
 *
 * Constructs event binning tool using command line arguments for user
 * parameter setting.
 ***************************************************************************/
ctbin::ctbin(int argc, char *argv[]) : 
       ctobservation(CTBIN_NAME, VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Event binning tool.
 *
 * Constructs event binning tool from another event binning tool.
 ***************************************************************************/
ctbin::ctbin(const ctbin& app) : ctobservation(app)
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
 *
 * Destructs event binning tool.
 ***************************************************************************/
ctbin::~ctbin(void)
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
 * @param[in] app Event binning tool.
 * @return Event binning tool.
 *
 * Assigns event binning tool.
 ***************************************************************************/
ctbin& ctbin::operator=(const ctbin& app)
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
 * @brief Clear event binning tool
 *
 * Clears event binning tool.
 ***************************************************************************/
void ctbin::clear(void)
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
 * @brief Run the event binning tool
 *
 * Gets the user parameters and loops over all CTA observations in the
 * observation container to bin the events into a single counts cube. All
 * observations in the observation container that do not contain CTA event
 * lists will be skipped.
 ***************************************************************************/
void ctbin::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Initialise sky direction cache for stacking
    init_sky_dir_cache();

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write header into logger
    log_header1(TERSE, gammalib::number("Find unbinned observation",
                m_obs.size()));

    // Find all unbinned CTA observations in m_obs
    std::vector<GCTAObservation*> obs_list(0);
    for (GCTAObservation* obs = first_unbinned_observation(); obs != NULL;
            obs = next_unbinned_observation()) {

        // Push observation into list
        obs_list.push_back(obs);

        // Write message
        std::string msg = " Including unbinned "+obs->instrument()+
                          " observation";
        log_string(NORMAL, msg);

    } // endfor: looped over all unbinned observations

    // Set number of relevant observations
    int nobs = obs_list.size();

    // Initialise arrays
    if (nobs > 0) {
        m_counts  = std::vector<GSkyMap>(nobs);
        m_weights = std::vector<GSkyMap>(nobs);
        m_cubes   = std::vector<GCTAEventCube>(nobs);
    }
    else {
        m_counts.clear();
        m_weights.clear();
        m_cubes.clear();
    }

    // Write header into logger
    log_header1(TERSE, gammalib::number("Bin observation", nobs));

    // Loop over all unbinned CTA observations in the container
    #pragma omp parallel for
    for (int i = 0; i < nobs; ++i) {

        // Get pointer to observation
        GCTAObservation* obs = obs_list[i];

        // Fill the cube
        GSkyMap counts = fill_cube(obs);

        // Set the counts cube weights
        GSkyMap weights = set_weights(obs);

        // Dispose events to free memory
        obs->dispose_events();

        // Post process cubes
        #pragma omp critical(ctbin_fill_cube)
        {
            // If stacking of cubes is requested then add up counts cube
            // and set weights in the first slot of the m_counts and
            // m_weights arrays
            if (m_stack) {
                if (i == 0) {
                    m_counts[0]  = counts;
                    m_weights[0] = weights;
                }
                else {
                    for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {
                        for (int pixel = 0; pixel < counts.npix(); ++pixel) {
                            m_counts[0](pixel, iebin) += counts(pixel, iebin);
                            if (weights(pixel, iebin) == 1.0) {
                                m_weights[0](pixel, iebin) = 1.0;
                            }
                        }
                    }
                }
            }

            // ... otherwise store counts and weights
        	else {
                m_counts[i]  = counts;
                m_weights[i] = weights;
            }

        } // end: omp critical section

    } // endfor: looped over observations

    // Combine GTIs of all observations and compute total ontime and livetime
    for (int i = 0; i < nobs; ++i) {

        // Get pointer to observation
        GCTAObservation* obs = obs_list[i];

        // Get Good Time Intervals for observation
        GGti gti = obs->gti();

        // If resulting GTI is not filled then copy over the reference from
        // the observation and use it for the resulting GTI
	    if (m_gti.is_empty()) {
            m_gti.reference(gti.reference());
	    }

	    // Append GTIs
	    m_gti.extend(gti);

	    // Update ontime and livetime
	    m_ontime   += obs->ontime();
	    m_livetime += obs->livetime();

    } // endfor: looped over all observations

    // Post process data for a stacked cube
    if (m_stack && nobs > 0) {

        // Compute livetime fraction per energy bin so that varying energy
        // thresholds are correctly taken into account in the weighting
        // array. Computational speed is not optimum, we have to see later
        // whether this is an issue.
        for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {
            double livetime_ebin = 0.0;
            for (int i = 0; i < obs_list.size(); ++i) {
                std::vector<bool> usage =
                     cube_layer_usage(m_ebounds, obs_list[i]->ebounds());
                if (usage[iebin]) {
                    livetime_ebin += obs_list[i]->livetime();
                }
            }
            if (m_livetime > 0.0) {
                double fraction = livetime_ebin / m_livetime;
                for (int pixel = 0; pixel < m_counts[0].npix(); ++pixel) {
                    if (m_weights[0](pixel, iebin) > 0.0) {
                        m_weights[0](pixel, iebin) *= fraction;
                    }
                }
            }
        }

        // Set a single cube in the observation container
        obs_cube_stacked();

    } // endif: stacked cube post processing

    // Post process data for non-stacked cubes
    else {

        // Loop over all relevant observations
        for (int i = 0; i < nobs; ++i) {

            // Get pointer to observation
            GCTAObservation* obs = obs_list[i];

            // Build counts cube
            m_cubes[i] = GCTAEventCube(m_counts[i], m_weights[i], m_ebounds,
                                       obs->gti());

            // Attach event cube in the correct slot of the observation
            // container
            obs->events(m_cubes[i]);

            // Construct unique filename for possible saving of the counts
            // cube into FITS files. Each observation has a unique pair of
            // instrument() and id() attributes, hence we construct a lower
            // case file name from these attributes, making sure that all
            // whitespaces are converted into "_" characters. Note that the
            // id() attribute is allowed to be blank, which can only happen
            // for a single observation.
            std::string inst = gammalib::tolower(obs->instrument());
            std::string id   = gammalib::tolower(obs->id());
            inst             = gammalib::replace_segment(inst, " ", "_");
            id               = gammalib::replace_segment(id, " ", "_");
            std::string filename = m_prefix + inst;
            if (!id.empty()) {
                filename += "_" + id;
            }
            filename += ".fits";

            // Set filename of observation
            obs->eventfile(filename);

        } // endfor: looped over relevant observations

    } // endelse: non-stacked cube post processing

    // Write resulting observation container into logger
    log_observations(NORMAL, m_obs, "Binned observation");

    // Optionally publish counts cube
    if (m_publish) {
        publish();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return event cube at index
 *
 * @param[in] index Cube index (0,...,cubes()-1).
 * @return Reference to event cube.
 *
 * @exception GException::out_of_range
 *            Event cube index out of range.
 *
 * Returns a reference to the event cube at the given index.
 ***************************************************************************/
const GCTAEventCube& ctbin::cube(const int& index) const
{
    // Compile option: raise an exception if index is out of range
    if (index < 0 || index >= cubes()) {
        throw GException::out_of_range(G_CUBE, "Cube index", index, cubes());
    }

    // Return cube
    return m_cubes[index];
}


/***********************************************************************//**
 * @brief Save counts cube
 *
 * Saves the counts cube. If the counts cube was stacked single FITS file
 * will be written. Otherwise each generated counts cube will be written
 * into a FITS file, and an observation definition XML file will be written
 * out that contains an observation for each counts cube.
 ***************************************************************************/
void ctbin::save(void)
{
    // Write header
    log_header1(TERSE, "Save counts cube");

    // If counts cube was stacked then save first observation
    if (m_stack) {

        // Save only if filename is valid and if there is at least one
        // observation
        if ((*this)["outobs"].is_valid() && m_obs.size() > 0) {

            // Get counts cube filename
            GFilename outobs = (*this)["outobs"].filename();

            // Get CTA observation from observation container
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

            // Save only if observation is valid
            if (obs != NULL) {

                // Log counts cube file name
                log_value(NORMAL, "Counts cube file", outobs.url());

                // Save cube
                obs->save(outobs, clobber());

            } // endif: observation was valid

        } // endif: outobs file was valid

    } // endif: counts cube was stacked

    // ... otherwise write out all relevant counts cubes and the observation
    // definition XML file
    else {

        // Loop over all observations
        for (int i = 0; i < m_obs.size(); ++i) {

            // Get CTA observation from observation container
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

            // Continue only if observation is valid
            if (obs != NULL) {

                // Retrieve filename
                std::string filename = obs->eventfile();

                // Save counts cube
                obs->save(filename, clobber());

            } // endif: observation was valid

        } // endfor: looped over all observations

        // Write XML file
        if ((*this)["outobs"].is_valid()) {

            // Get observation definition XML filename
            std::string outobs = (*this)["outobs"].filename();

            // Issue warning if output filename has no .xml suffix
            log_string(TERSE, warn_xml_suffix(outobs));

            // Save XML file
            m_obs.save(outobs);
        }

    } // endelse: counts cube were not stacked

    // Return
    return;
}


/***********************************************************************//**
 * @brief Publish counts cube
 *
 * @param[in] name Counts cube name.
 ***************************************************************************/
void ctbin::publish(const std::string& name)
{
    // Write header into logger
    log_header1(TERSE, "Publish counts cube");

    // Set default name is user name is empty
    std::string user_name(name);
    if (user_name.empty()) {
        user_name = CTBIN_NAME;
    }

    // Write counts cube name into logger
    log_value(NORMAL, "Counts cube name", user_name);

    // Publish first counts cube if it exists
    if (cubes() > 0) {
        m_counts[0].publish(user_name);
    }

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
void ctbin::init_members(void)
{
    // Initialise members
    m_usepnt  = false;
    m_stack   = true;
    m_prefix.clear();
    m_publish = false;
    m_chatter = static_cast<GChatter>(2);

    // Initialise protected members
    m_cubes.clear();
    m_counts.clear();
    m_weights.clear();
    m_ebounds.clear();
    m_gti.clear();
    m_ontime   = 0.0;
    m_livetime = 0.0;

    // Initialise cache members
    m_dirs.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctbin::copy_members(const ctbin& app)
{
    // Copy attributes
    m_usepnt  = app.m_usepnt;
    m_stack   = app.m_stack;
    m_prefix  = app.m_prefix;
    m_publish = app.m_publish;
    m_chatter = app.m_chatter;

    // Copy protected members
    m_cubes    = app.m_cubes;
    m_counts   = app.m_counts;
    m_weights  = app.m_weights;
    m_ebounds  = app.m_ebounds;
    m_gti      = app.m_gti;
    m_ontime   = app.m_ontime;
    m_livetime = app.m_livetime;

    // Copy cache members
    m_dirs = app.m_dirs;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctbin::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. Most parameters are only required if no observation exists so
 * far in the observation container. In this case, a single CTA observation
 * will be added to the container, using the definition provided in the
 * parameter file.
 ***************************************************************************/
void ctbin::get_parameters(void)
{
    // Setup observations from "inobs" parameter. Do not request response
    // information and do not accept counts cubes.
    setup_observations(m_obs, false, true, false);

    // Create an event cube to query task parameters
    GCTAEventCube cube = this->ctool::create_cube(m_obs);

    // Get energy boundaries
    m_ebounds = cube.ebounds();

    // Get stack and prefix parameters
    m_stack  = (*this)["stack"].boolean();
    m_prefix = (*this)["prefix"].string();

    // Get remaining parameters
    m_publish = (*this)["publish"].boolean();
    m_chatter = static_cast<GChatter>((*this)["chatter"].integer());

    // If needed later, query output filename now
    if (read_ahead()) {
        (*this)["outobs"].query();
    }

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
 * @brief Initialise sky direction cache for cube stack
 *
 * Initialises a cache that holds the sky directions of each counts cube
 * pixel for a stacked analysis. This will speed up the computations later.
 *
 * If no stacking is requested, the method does nothing.
 ***************************************************************************/
void ctbin::init_sky_dir_cache(void)
{
    // Continue only is stacking is requested
    if (m_stack) {

        // Create a sky map from task parameters
        GSkyMap map = this->ctool::create_map(m_obs);

        // Initialise sky direction
        m_dirs.resize(map.npix());

        // Compute sky direction for each pixel
        for (int pixel = 0; pixel < map.npix(); ++pixel) {
            m_dirs[pixel] = map.inx2dir(pixel);
        }

    } // endif: stacking was requested

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create counts cube sky map for current observation
 *
 * @param[in] obs CTA observation.
 *
 * Creates a counts cube sky map for current observations.
 ***************************************************************************/
GSkyMap ctbin::create_cube(const GCTAObservation* obs)
{
    // Read coordinate system and projection
    std::string coordsys = (*this)["coordsys"].string();
    std::string proj     = (*this)["proj"].string();

    // Read sky map reference
    double xref   = 0.0;
    double yref   = 0.0;
    bool   usepnt = (*this)["usepnt"].boolean();
    if (!usepnt) {
        xref = (*this)["xref"].real();
        yref = (*this)["yref"].real();
    }

    // Read sky map pixel size and number of bins
    double binsz = (*this)["binsz"].real();
    int    nxpix = (*this)["nxpix"].integer();
    int    nypix = (*this)["nypix"].integer();

    // Read energy boundaries
    GEbounds ebounds = create_ebounds();

    // If requested, get pointing from observations
    if (usepnt) {

        // Get pointing direction for observation
        GSkyDir pnt = obs->pointing().dir();

        // Set xref/yref based on the coordinate system
        if (gammalib::toupper(coordsys) == "GAL") {
            xref = pnt.l_deg();
            yref = pnt.b_deg();
        }
        else {
            xref = pnt.ra_deg();
            yref = pnt.dec_deg();
        }

    } // endif: got pointing from observation

    // Initialise counts cube
    GSkyMap cube = GSkyMap(proj, coordsys, xref, yref, -binsz, binsz,
                           nxpix, nypix, ebounds.size());


    // Return cube
    return cube;
}


/***********************************************************************//**
 * @brief Fill events into counts cube
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No event list or valid RoI found in observation.
 *
 * Fills the events from an event list into the counts cube.
 ***************************************************************************/
GSkyMap ctbin::fill_cube(const GCTAObservation* obs)
{
    // Make sure that the observation holds a CTA event list. If this is
    // not the case then throw an exception.
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>
                                  (obs->events());
    if (events == NULL) {
        std::string msg = "CTA Observation does not contain an event "
                          "list. An event list is needed to fill the "
                          "counts cube.";
        throw GException::invalid_value(G_FILL_CUBE, msg);
    }

    // Get the RoI
    const GCTARoi& roi = obs->roi();

    // Check for RoI sanity
    if (!roi.is_valid()) {
        std::string msg = "No valid RoI found in input observation "
                          "\""+obs->name()+"\". Run ctselect to specify "
                          "an RoI for this observation before running "
                          "ctbin.";
        throw GException::invalid_value(G_FILL_CUBE, msg);
    }

    // Get counts cube usage flags
    std::vector<bool> usage = cube_layer_usage(m_ebounds, obs->ebounds());

    // Initialise binning statistics
    int num_outside_roi  = 0;
    int num_invalid_wcs  = 0;
    int num_outside_map  = 0;
    int num_outside_ebds = 0;
    int num_in_map       = 0;

    // Create counts cube
    GSkyMap counts = create_cube(obs);

    // Extract counts cube boundaries
    double pixel_x_max = counts.nx() - 0.5;
    double pixel_y_max = counts.ny() - 0.5;

    // Fill counts sky map
    for (int i = 0; i < events->size(); ++i) {

        // Get event
        const GCTAEventAtom* event = (*events)[i];

        // Determine event sky direction
        GCTAInstDir* inst = (GCTAInstDir*)&(event->dir());
        GSkyDir      dir  = inst->dir();

        // Skip event if it is outside the RoI
        if (roi.centre().dir().dist_deg(dir) > roi.radius()) {
            num_outside_roi++;
            continue;
        }

        // Determine sky pixel
        GSkyPixel pixel;
        try {
            pixel = counts.dir2pix(dir);
        }
        catch (std::exception &e) {
            num_invalid_wcs++;
            continue;
        }

        // Skip event if corresponding counts cube pixel is outside the
        // counts cube map range
        if (pixel.x() < -0.5 || pixel.x() > pixel_x_max ||
            pixel.y() < -0.5 || pixel.y() > pixel_y_max) {
            num_outside_map++;
            continue;
        }

        // Determine counts cube energy bin
        int iebin = m_ebounds.index(event->energy());

        // Skip event if the corresponding counts cube energy bin is not
        // fully contained in the event list energy range. This avoids
        // having partially filled bins.
        if (iebin == -1 || !usage[iebin]) {
            num_outside_ebds++;
            continue;
        }

        // Fill event in skymap
        counts(pixel, iebin) += 1.0;

        // Increment number of maps
        num_in_map++;

    } // endfor: looped over all events

    // Log filling results
    #pragma omp critical(ctbin_fill_cube)
    {
        log_header3(TERSE, get_obs_header(obs));
        log_value(NORMAL, "Events in list", obs->events()->size());
        log_value(NORMAL, "Events in cube", num_in_map);
        log_value(NORMAL, "Events outside RoI", num_outside_roi);
        log_value(NORMAL, "Events with invalid WCS", num_invalid_wcs);
        log_value(NORMAL, "Events outside cube area", num_outside_map);
        log_value(NORMAL, "Events outside energy bins", num_outside_ebds);
    }

    // Return
    return counts;
}


/***********************************************************************//**
 * @brief Set counts cube weights for a given observation
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No event list or valid RoI found in observation.
 *
 * Sets the counts cube weights for all bins that are considered for the
 * specific observation to unity.
 ***************************************************************************/
GSkyMap ctbin::set_weights(const GCTAObservation* obs)
{
    // Get the RoI
    const GCTARoi& roi = obs->roi();

    // Check for RoI sanity
    if (!roi.is_valid()) {
        std::string msg = "No valid RoI found in input observation "
                          "\""+obs->name()+"\". Run ctselect to specify "
                          "an RoI for this observation before running "
                          "ctbin.";
        throw GException::invalid_value(G_SET_WEIGHTS, msg);
    }

    // Get the RoI centre and radius in radians
    GSkyDir roi_centre = roi.centre().dir();
    double  roi_radius = roi.radius() * gammalib::deg2rad;

    // Get counts cube usage flags
    std::vector<bool> usage = cube_layer_usage(m_ebounds, obs->ebounds());

    // Create weights cube
    GSkyMap weights = create_cube(obs);

    // Set sky direction cache availability flag
    bool cache = (m_stack && m_dirs.size() == weights.npix());

    // Loop over all pixels in counts cube
    for (int pixel = 0; pixel < weights.npix(); ++pixel) {

        // Skip pixel if it is outside the RoI
        if (cache) {
            if (m_dirs[pixel].dist(roi_centre) > roi_radius) {
                continue;
            }
        }
        else {
            if (roi_centre.dist(weights.inx2dir(pixel)) > roi_radius) {
                continue;
            }
        }

        // Loop over all energy layers of counts cube
        for (int iebin = 0; iebin < m_ebounds.size(); ++iebin){

            // Skip energy layer if the usage flag is false
            if (!usage[iebin]) {
                continue;
            }

            // Signal that bin was filled
            weights(pixel, iebin) = 1.0;

        } // endfor: looped over energy layers of counts cube

    } // endfor: looped over pixels of counts cube

    // Return
    return weights;
}


/***********************************************************************//**
 * @brief Create output observation container.
 *
 * Creates an output observation container that combines all input CTA
 * observation into a single stacked observation. All non-CTA observations
 * and all binned CTA observations that were present in the observation
 * container are append to the observation container so that they can
 * be used by other tools. The method furthermore conserves any response
 * information in case that a single CTA observation is provided to support
 * binned analysis.
 ***************************************************************************/
void ctbin::obs_cube_stacked(void)
{
    // Set event cube in first slot
    m_cubes[0] = GCTAEventCube(m_counts[0], m_weights[0], m_ebounds, m_gti);

    // If we have only a single CTA observation in the container, then
    // keep that observation and just attach the event cube to it. Reset
    // the filename, otherwise we still will have the old event filename
    // in the log file.
    if (m_obs.size() == 1) {

        // Attach event cube to CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);
        if (obs != NULL) {

            // Change the event type if we had an unbinned observation
            if (obs->eventtype() == "EventList") {

                // Assign cube to the observation
                obs->events(m_cubes[0]);

            }

            // ... otherwise the input observation was binned and hence
            // skipped. In that case we simply append an empty counts
            // cube
            else {

                // Create empty counts cube
                GCTAEventCube cube(m_counts[0], m_weights[0], m_ebounds,
                                   obs->gti());

                // Assign empty cube to observation
                obs->events(cube);

            }

            // Reset file name
            obs->eventfile("");

        } // endif: obervation was valid

    } // endif: we only had one observation in the container

    // ... otherwise put a single CTA observation in container and
    // append all observations that have not been used in the binning
    // to the container
    else {

        // Allocate observation container
        GObservations container;

        // Allocate CTA observation.
        GCTAObservation obs;

        // Attach event cube to CTA observation
        obs.events(m_cubes[0]);

        // Compute average pointing direction for all CTA event lists
        double ra     = 0.0;
        double dec    = 0.0;
        double number = 0.0;
        for (int i = 0; i < m_obs.size(); ++i) {
            GCTAObservation* cta = dynamic_cast<GCTAObservation*>(m_obs[i]);
            if ((cta != NULL) && (cta->eventtype() == "EventList")) {
                ra     += cta->pointing().dir().ra();
                dec    += cta->pointing().dir().dec();
                number += 1.0;
            }
        }
        if (number > 0.0) {
            ra  /= number;
            dec /= number;
        }
        GSkyDir dir;
        dir.radec(ra, dec);
        GCTAPointing pointing(dir);

        // Compute deadtime correction
        double deadc = (m_ontime > 0.0) ? m_livetime / m_ontime : 0.0;

        // Set CTA observation attributes
        obs.pointing(pointing);
        obs.ra_obj(dir.ra_deg());   //!< Dummy
        obs.dec_obj(dir.dec_deg()); //!< Dummy
        obs.ontime(m_ontime);
        obs.livetime(m_livetime);
        obs.deadc(deadc);

        // Set models in observation container
        container.models(m_obs.models());

        // Append CTA observation
        container.append(obs);

        // Copy over all remaining non-CTA observations
        for (int i = 0; i < m_obs.size(); ++i) {
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
            if (obs == NULL) {
                container.append(*m_obs[i]);
            }
            else if (obs->eventtype() != "EventList") {
                container.append(*m_obs[i]);
            }
        }

        // Set observation container
        m_obs = container;

    } // endelse: there was not a single CTA observation

    // Return
    return;
}
