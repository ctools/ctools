/***************************************************************************
 *                        ctbin - Event binning tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
#define G_GET_PARAMETERS          "ctbin::get_parameters()"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
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
 * This method creates an instance of the class by copying an existing
 * observations container.
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
 * @brief Clear instance
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
 * @brief Bin the event data
 *
 * This method loops over all observations found in the observation conatiner
 * and bins all events from the event list(s) into counts map(s). Note that
 * each event list is binned in a separate counts map, hence no summing of
 * events is done.
 ***************************************************************************/
void ctbin::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Observations");
        }
        else {
            log.header1("Observation");
        }
        log << m_obs << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Bin observations");
        }
        else {
            log.header1("Bin observation");
        }
    }

    // Loop over all observations in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Continue only if observation is a CTA observation
        if (obs != NULL) {

            // Write header for observation
            if (logTerse()) {
                if (obs->name().length() > 1) {
                    log.header3("Observation "+obs->name());
                }
                else {
                    log.header3("Observation");
                }
            }

            // Fill the cube
            fill_cube(obs);

            // Dispose events to free memory
            obs->dispose_events();

        } // endif: CTA observation found

    } // endfor: looped over observations

    // Set a single cube in the observation container
    obs_cube();

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Binned observations");
        }
        else {
            log.header1("Binned observation");
        }
        log << m_obs << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save counts cube
 *
 * This method saves the counts cube into a FITS file.
 ***************************************************************************/
void ctbin::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Save observations");
        }
        else {
            log.header1("Save observation");
        }
    }

    // Get output filename
    m_outcube = (*this)["outcube"].filename();

    // Get CTA observation from observation container
    GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

    // Save only if observation is valid
    if (obs != NULL) {
        obs->save(m_outcube, clobber());
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
    m_outcube.clear();
    m_usepnt   = false;

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
    m_outcube  = app.m_outcube;
    m_usepnt   = app.m_usepnt;

    // Copy protected members
    m_obs        = app.m_obs;
    m_cube       = app.m_cube;
    m_ebounds    = app.m_ebounds;
    m_gti        = app.m_gti;
    m_ontime     = app.m_ontime;
    m_livetime   = app.m_livetime;

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
 *  @exception GException::invalid_value
 *            Parameter "inobs" is required for ctbin.
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. Most parameters are only required if no observation exists so
 * far in the observation container. In this case, a single CTA observation
 * will be added to the container, using the definition provided in the
 * parameter file.
 ***************************************************************************/
void ctbin::get_parameters(void)
{
    // If there are no observations in container then load them via user parameters
    if (m_obs.size() == 0) {

        // Throw exception if no infile is given, since this tool needs an observation
        // including events
        if ((*this)["inobs"].filename()=="NONE" || (*this)["inobs"].filename() == "") {

            std::string msg = "Parameter \"inobs\" is required to be given in ctbin."
                            "Specify a vaild observation definition (XML or FITS) file to proceed";
            throw GException::invalid_value(G_GET_PARAMETERS, msg);
        }

        // Build observation container without response (not needed)
        m_obs = get_observations(false);

    } // endif: there was no observation in the container

    // Get parameters
    m_usepnt = (*this)["usepnt"].boolean();

    // Initialise event cube
    GCTAEventCube cube;

    // Check if pointing should be used
    if (m_usepnt) {

        // buld cube from user parameters and pointing information
        cube = get_cube(m_obs);

    } // endif: pointing used as map center

    else {

        // Build cube from user parameters
        cube = get_cube();

    } // endelse: m_usepnt was false

    // Get the skymap
    m_cube = cube.map();

    // Get energy boundaries
    m_ebounds  = cube.ebounds();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outcube = (*this)["outcube"].filename();
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
                              "list. Event list information is needed to "
                              "fill the counts map.";
            throw GException::invalid_value(G_FILL_CUBE, msg);
        }

        // Get the RoI
        const GCTARoi& roi = events->roi();

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
            if (index == -1) {
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
 * observation into a single cube-style observation. All non CTA observations
 * present in the observation container are kept. The method furthermore
 * conserves any response information in case that a single CTA observation
 * is provided. This supports the original binned analysis.
 ***************************************************************************/
void ctbin::obs_cube(void)
{
    // If we have only a single CTA observation in the container, then
    // keep that observation and just attach the event cube to it
    if (m_obs.size() == 1) {

        // Attach event cube to CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);
        if (obs != NULL) {
            obs->events(this->cube());
        }

    }

    // ... otherwise put a single CTA observation in container
    else {

        // Allocate observation container
        GObservations container;

        // Allocate CTA observation.
        GCTAObservation obs;

        // Attach event cube to CTA observation
        obs.events(this->cube());

        // Set map centre as pointing
        GSkyPixel    pixel(0.5*double(m_cube.nx()), 0.5*double(m_cube.ny()));
        GSkyDir      centre = m_cube.pix2dir(pixel);
        GCTAPointing pointing(centre);

        // Compute deadtime correction
        double deadc = (m_ontime > 0.0) ? m_livetime / m_ontime : 0.0;

        // Set CTA observation attributes
        obs.pointing(pointing);
        obs.obs_id(0);
        obs.ra_obj(centre.ra_deg());   //!< Dummy
        obs.dec_obj(centre.dec_deg()); //!< Dummy
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
        }

        // Set observation container
        m_obs = container;

    } // endelse: there was not a single CTA observation

    // Return
    return;
}
