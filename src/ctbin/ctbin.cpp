/***************************************************************************
 *                        ctbin - Event binning tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2017 by Juergen Knoedlseder                         *
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

/* __ Method name definitions ____________________________________________ */
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
ctbin::ctbin(void) : ctobservation(CTBIN_NAME, CTBIN_VERSION)
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
       ctobservation(CTBIN_NAME, CTBIN_VERSION, obs)
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
       ctobservation(CTBIN_NAME, CTBIN_VERSION, argc, argv)
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

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write header into logger
    log_header1(TERSE, gammalib::number("Bin observation", m_obs.size()));

    // Loop over all unbinned CTA observations in the container
    for (GCTAObservation* obs = first_unbinned_observation(); obs != NULL;
         obs = next_unbinned_observation()) {

        // Fill the cube
        fill_cube(obs);

        // Set the counts cube weights
        set_weights(obs);

        // Dispose events to free memory
        obs->dispose_events();

    } // endfor: looped over observations

    // Build event cube (needs to come before obs_cube() since this method
    // relies on correct setting of m_cube)
    m_cube = GCTAEventCube(m_counts, m_weights, m_ebounds, m_gti);

    // Set a single cube in the observation container
    obs_cube();

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
 * @brief Save counts cube
 *
 * Saves the counts cube into a FITS file.
 ***************************************************************************/
void ctbin::save(void)
{
    // Write header
    log_header1(TERSE, "Save counts cube");

    // Get counts cube filename
    m_outcube = (*this)["outcube"].filename();

    // Save only if filename is not empty and if there is at least one
    // observation
    if (!m_outcube.is_empty() && m_obs.size() > 0) {

        // Get CTA observation from observation container
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

        // Save only if observation is valid
        if (obs != NULL) {
        
            // Log counts cube file name
            log_value(NORMAL, "Counts cube file", m_outcube.url());
            
            // Save cube
            obs->save(m_outcube, clobber());

        } // endif: observation was valid

    } // endif: outcube file was valid

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

    // Publish counts cube
    m_counts.publish(user_name);

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
    m_outcube.clear();
    m_usepnt  = false;
    m_publish = false;
    m_chatter = static_cast<GChatter>(2);

    // Initialise protected members
    m_counts.clear();
    m_weights.clear();
    m_ebounds.clear();
    m_gti.clear();
    m_cube.clear();
    m_ontime   = 0.0;
    m_livetime = 0.0;

    // Set CTA reference time
    m_gti.reference(m_cta_ref);

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
    m_outcube = app.m_outcube;
    m_usepnt  = app.m_usepnt;
    m_publish = app.m_publish;
    m_chatter = app.m_chatter;

    // Copy protected members
    m_counts    = app.m_counts;
    m_weights   = app.m_weights;
    m_ebounds   = app.m_ebounds;
    m_gti       = app.m_gti;
    m_cube      = app.m_cube;
    m_ontime    = app.m_ontime;
    m_livetime  = app.m_livetime;

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

    // Create an event cube based on task parameters
    GCTAEventCube cube = create_cube(m_obs);

    // Get the skymap from the cube and initialise all counts cube bins and
    // weights to zero
    m_counts  = cube.counts();
    m_counts  = 0.0;
    m_weights = m_counts;

    // Get energy boundaries
    m_ebounds = cube.ebounds();

    // Get remaining parameters
    m_publish = (*this)["publish"].boolean();
    m_chatter = static_cast<GChatter>((*this)["chatter"].integer());

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outcube = (*this)["outcube"].filename();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
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
void ctbin::fill_cube(GCTAObservation* obs)
{
    // Make sure that the observation holds a CTA event list. If this
    // is not the case then throw an exception.
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>
                                  (obs->events());
    if (events == NULL) {
        std::string msg = "CTA Observation does not contain an event "
                          "list. An event list is needed to fill the "
                          "counts cube.";
        throw GException::invalid_value(G_FILL_CUBE, msg);
    }

    // Get the RoI
    const GCTARoi& roi = events->roi();

    // Check for RoI sanity
    if (!roi.is_valid()) {
        std::string msg = "No valid RoI found in input observation "
                          "\""+obs->name()+"\". Run ctselect to specify "
                          "an RoI for this observation before running "
                          "ctbin.";
        throw GException::invalid_value(G_FILL_CUBE, msg);
    }

    // Get counts cube usage flags
    std::vector<bool> usage = cube_layer_usage(m_ebounds, events->ebounds());

    // Initialise binning statistics
    int num_outside_roi  = 0;
    int num_outside_map  = 0;
    int num_outside_ebds = 0;
    int num_in_map       = 0;

    // Fill counts sky map
    for (int i = 0; i < events->size(); ++i) {

        // Get event
        const GCTAEventAtom* event = (*events)[i];

        // Determine sky pixel
        GCTAInstDir* inst  = (GCTAInstDir*)&(event->dir());
        GSkyDir      dir   = inst->dir();
        GSkyPixel    pixel = m_counts.dir2pix(dir);

        // Skip event if corresponding counts cube pixel is outside RoI
        if (roi.centre().dir().dist_deg(dir) > roi.radius()) {
            num_outside_roi++;
            continue;
        }

        // Skip event if corresponding counts cube pixel is outside the
        // counts cube map range
        if (pixel.x() < -0.5 || pixel.x() > (m_counts.nx()-0.5) ||
            pixel.y() < -0.5 || pixel.y() > (m_counts.ny()-0.5)) {
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
        m_counts(pixel, iebin) += 1.0;
        num_in_map++;

    } // endfor: looped over all events

    // Append GTIs
    m_gti.extend(events->gti());

    // Update ontime and livetime
    m_ontime   += obs->ontime();
    m_livetime += obs->livetime();

    // Log filling results
    log_value(NORMAL, "Events in list", obs->events()->size());
    log_value(NORMAL, "Events in cube", num_in_map);
    log_value(NORMAL, "Event bins outside RoI", num_outside_roi);
    log_value(NORMAL, "Events outside cube area", num_outside_map);
    log_value(NORMAL, "Events outside energy bins", num_outside_ebds);

    // Return
    return;
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
void ctbin::set_weights(GCTAObservation* obs)
{
    // Make sure that the observation holds a CTA event list. If this
    // is not the case then throw an exception.
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>
                                  (obs->events());
    if (events == NULL) {
        std::string msg = "CTA Observation does not contain an event "
                          "list. An event list is needed to fill the "
                          "counts cube.";
        throw GException::invalid_value(G_SET_WEIGHTS, msg);
    }

    // Get the RoI
    const GCTARoi& roi = events->roi();

    // Check for RoI sanity
    if (!roi.is_valid()) {
        std::string msg = "No valid RoI found in input observation "
                          "\""+obs->name()+"\". Run ctselect to specify "
                          "an RoI for this observation before running "
                          "ctbin.";
        throw GException::invalid_value(G_SET_WEIGHTS, msg);
    }

    // Get counts cube usage flags
    std::vector<bool> usage = cube_layer_usage(m_ebounds, events->ebounds());

    // Loop over all pixels in counts cube
    for (int pixel = 0; pixel < m_counts.npix(); ++pixel) {

        // Get pixel sky direction
        GSkyDir dir = m_counts.inx2dir(pixel);

        // Skip pixel if it is outside the RoI
        if (roi.centre().dir().dist_deg(dir) > roi.radius()) {
            continue;
        }

        // Loop over all energy layers of counts cube
        for (int iebin = 0; iebin < m_ebounds.size(); ++iebin){

            // Skip energy layer if the usage flag is false
            if (!usage[iebin]) {
                continue;
            }
            // Signal that bin was filled
            m_weights(pixel, iebin) = 1.0;

        } // endfor: looped over energy layers of counts cube

    } // endfor: looped over pixels of counts cube

    // Return
    return;
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
void ctbin::obs_cube(void)
{
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
                obs->events(this->cube());

            }

            // ... otherwise the input observation was binned and hence
            // skipped. In that case we simply append an empty counts
            // cube
            else {

                // Create empty counts cube
                GCTAEventCube cube(m_counts, m_weights, m_ebounds, obs->gti());

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
        obs.events(this->cube());

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
        obs.obs_id(0);
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
