/***************************************************************************
 *                        ctbin - Event binning tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
#define G_FILL_CUBE                      "ctbin::fill_cube(GCTAObservation*)"
#define G_GET_PARAMETERS                            "ctbin::get_parameters()"

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
 *
 * Constructs an empty ctbin tool.
 ***************************************************************************/
ctbin::ctbin(void) : ctool(CTBIN_NAME, CTBIN_VERSION)
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
 * Constructs a ctbin tool from an observation container.
 ***************************************************************************/
ctbin::ctbin(const GObservations& obs) : ctool(CTBIN_NAME, CTBIN_VERSION)
{
    // Initialise members
    init_members();

    // Set observations
    m_obs = obs;

    // Return
    return;
}



/***********************************************************************//**
 * @brief Command line constructor
 *
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 *
 * Constructs ctbin tool from command line arguments.
 ***************************************************************************/
ctbin::ctbin(int argc, char *argv[]) : 
       ctool(CTBIN_NAME, CTBIN_VERSION, argc, argv)
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
 *
 * Constructs ctbin tool from another ctbin instance.
 ***************************************************************************/
ctbin::ctbin(const ctbin& app) : ctool(app)
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
 * Destructs ctbin tool.
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
 * @param[in] app Application.
 * @return ctbin tool.
 *
 * Assigns ctbin tool.
 ***************************************************************************/
ctbin& ctbin::operator=(const ctbin& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->ctool::operator=(app);

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
 * @brief Clear ctbin tool
 *
 * Clears ctbin tool.
 ***************************************************************************/
void ctbin::clear(void)
{
    // Free members
    free_members();
    this->ctool::free_members();
    this->GApplication::free_members();

    // Initialise members
    this->GApplication::init_members();
    this->ctool::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Run the ctbin tool
 *
 * Gets the user parameters and loops over all CTA observations in the
 * observation container to bin the events into a single counts cube.
 ***************************************************************************/
void ctbin::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        log.header1(gammalib::number("Observation", m_obs.size()));
        log << m_obs << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1(gammalib::number("Bin observation", m_obs.size()));
    }

    // Loop over all observations in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Write header for observation
        if (logTerse()) {
            std::string header = m_obs[i]->instrument() + " observation";
            if (!m_obs[i]->name().empty()) {
                header += " \"" + m_obs[i]->name() + "\"";
            }
            if (!m_obs[i]->id().empty()) {
                header += " (id=" + m_obs[i]->id() +")";
            }
            log.header3(header);
        }

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Skip observation if it's not CTA
        if (obs == NULL) {
            if (logTerse()) {
                log << " Skipping ";
                log << m_obs[i]->instrument();
                log << " observation" << std::endl;
            }
            continue;
        }

        // Skip observation if we have a binned observation
        if (obs->eventtype() == "CountsCube") {
            if (logTerse()) {
                log << " Skipping binned ";
                log << obs->instrument();
                log << " observation" << std::endl;
            }
            continue;
        }

        // Fill the cube
        fill_cube(obs);

        // Dispose events to free memory
        obs->dispose_events();

    } // endfor: looped over observations

    // Set a single cube in the observation container
    obs_cube();

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        log.header1(gammalib::number("Binned observation", m_obs.size()));
        log << m_obs << std::endl;
    }

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
    if (logTerse()) {
        log << std::endl;
        log.header1("Save counts cube");
    }

    // Get counts cube filename
    m_outcube = (*this)["outcube"].filename();

    // Save only if filename is non-empty
    if (!m_outcube.is_empty()) {

        // Get CTA observation from observation container
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

        // Save only if observation is valid
        if (obs != NULL) {
        
            // Log filename
            if (logTerse()) {
                log << "Save counts cube into file \""+m_outcube+"\".";
                log << std::endl;
            }
            
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
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Publish counts cube");
    }

    // Set default name is user name is empty
    std::string user_name(name);
    if (user_name.empty()) {
        user_name = CTBIN_NAME;
    }

    // Log filename
    if (logTerse()) {
        log << "Publish \""+user_name+"\" counts cube." << std::endl;
    }

    // Publish counts cube
    m_cube.publish(user_name);

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

    // Initialise protected members
    m_obs.clear();
    m_cube.clear();
    m_ebounds.clear();
    m_gti.clear();
    m_ontime   = 0.0;
    m_livetime = 0.0;

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

    // Copy protected members
    m_obs      = app.m_obs;
    m_cube     = app.m_cube;
    m_ebounds  = app.m_ebounds;
    m_gti      = app.m_gti;
    m_ontime   = app.m_ontime;
    m_livetime = app.m_livetime;

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
    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {

        // Throw exception if no input observation file is given
        require_inobs(G_GET_PARAMETERS);

        // Throw exception if counts cube is given
        require_inobs_nocube(G_GET_PARAMETERS);

        // Get observation container without response (not needed)
        m_obs = get_observations(false);

    } // endif: there was no observation in the container

    // Create an event cube based on task parameters
    GCTAEventCube cube = create_cube(m_obs);

    // Get the skymap from the cube and initialise all pixels to zero
    m_cube = cube.map();
    m_cube = 0.0;

    // Get energy boundaries
    m_ebounds  = cube.ebounds();

    // Get remaining parameters
    m_publish = (*this)["publish"].boolean();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outcube = (*this)["outcube"].filename();
    }

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill events into counts cube
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No event list found in observation.
 *
 * Fills the events from an event list in the counts cube setup by init_cube.
 ***************************************************************************/
void ctbin::fill_cube(GCTAObservation* obs)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Make sure that the observation holds a CTA event list. If this
        // is not the case then throw an exception.
        const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(obs->events());
        if (events == NULL) {
            std::string msg = "CTA Observation does not contain an event "
                              "list. An event list is needed to fill the "
                              "counts cube.";
            throw GException::invalid_value(G_FILL_CUBE, msg);
        }

        // Get the RoI
        const GCTARoi& roi = events->roi();

        // Get the ebounds
        const GEbounds& obs_ebounds = events->ebounds();

        // Check for RoI sanity
        if (!roi.is_valid()) {
            std::string msg = "No valid RoI found in input observation "
                              "\""+obs->name()+"\". Run ctselect to specify "
                              "an RoI for this observation before running "
                              "ctbin.";
            throw GException::invalid_value(G_FILL_CUBE, msg);
        }

        // Initialise binning statistics
        int num_outside_roi  = 0;
        int num_outside_map  = 0;
        int num_outside_ebds = 0;
        int num_in_map       = 0;

        // Fill sky map
        for (int i = 0; i < events->size(); ++i) {

            // Get event
            const GCTAEventAtom* event = (*events)[i];

            // Determine sky pixel
            GCTAInstDir* inst  = (GCTAInstDir*)&(event->dir());
            GSkyDir      dir   = inst->dir();
            GSkyPixel    pixel = m_cube.dir2pix(dir);

            // Skip if pixel is outside RoI
            if (roi.centre().dir().dist_deg(dir) > roi.radius()) {
                num_outside_roi++;
                continue;
            }

            // Skip if pixel is out of range
            if (pixel.x() < -0.5 || pixel.x() > (m_cube.nx()-0.5) ||
                pixel.y() < -0.5 || pixel.y() > (m_cube.ny()-0.5)) {
                num_outside_map++;
                continue;
            }

            // Determine energy bin. Skip if we are outside the energy range
            int index = m_ebounds.index(event->energy());
            if (index == -1 ||
                !obs_ebounds.contains(m_ebounds.emin(index)+g_energy_margin,
                                      m_ebounds.emax(index)-g_energy_margin)) {
                num_outside_ebds++;
                continue;
            }

            // Fill event in skymap
            m_cube(pixel, index) += 1.0;
            num_in_map++;

        } // endfor: looped over all events

        // Append GTIs
        m_gti.extend(events->gti());

        // Update ontime and livetime
        m_ontime   += obs->ontime();
        m_livetime += obs->livetime();

        // Log filling results
        if (logTerse()) {
            log << gammalib::parformat("Events in list");
            log << obs->events()->size() << std::endl;
            log << gammalib::parformat("Events in cube");
            log << num_in_map << std::endl;
            log << gammalib::parformat("Event bins outside RoI");
            log << num_outside_roi << std::endl;
            log << gammalib::parformat("Events outside cube area");
            log << num_outside_map << std::endl;
            log << gammalib::parformat("Events outside energy bins");
            log << num_outside_ebds << std::endl;
        }

        // Log cube
        if (logExplicit()) {
            log.header1("Counts cube");
            log << m_cube << std::endl;
        }

    } // endif: observation was valid

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
                GCTAEventCube cube(m_cube, m_ebounds, obs->gti());

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
