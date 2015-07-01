/***************************************************************************
 *                       ctskymap - Sky mapping tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2015 by Juergen Knoedlseder                         *
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
 * @file ctskymap.cpp
 * @brief Sky mapping tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctskymap.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_INIT_MAP                 "ctskymap::init_map(GCTAObservation* obs)"
#define G_BIN_EVENTS                 "ctskymap::bin_events(GCTAObservation*)"
#define G_GET_PARAMETERS                         "ctskymap::get_parameters()"

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
ctskymap::ctskymap(void) : ctool(CTSKYMAP_NAME, CTSKYMAP_VERSION)
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
ctskymap::ctskymap(const GObservations& obs) :
          ctool(CTSKYMAP_NAME, CTSKYMAP_VERSION)
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
ctskymap::ctskymap(int argc, char *argv[]) : 
          ctool(CTSKYMAP_NAME, CTSKYMAP_VERSION, argc, argv)
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
ctskymap::ctskymap(const ctskymap& app) : ctool(app)
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
ctskymap::~ctskymap(void)
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
ctskymap& ctskymap::operator=(const ctskymap& app)
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
void ctskymap::clear(void)
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
 * @brief Creates sky maps from data
 *
 * This method is the main code. It
 * (1) reads task parameters from the par file
 * (2) initialises the sky maps
 * (3) loops over all observations to add its events to the sky map
 ***************************************************************************/
void ctskymap::run(void)
{
    // Initialise statistics
    int num_obs = 0;

    // Switch screen logging on in debug mode
    if (logDebug()) {
        log.cout(true);
    }

    // Get parameters
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
            log.header1("Map observations");
        }
        else {
            log.header1("Map observation");
        }
    }

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Skip observation if it's not CTA
        if (obs == NULL) {
            if (logTerse()) {
                log << "Warning: Skipping "+m_obs[i]->instrument();
                log << " observation \"";
                log << m_obs[i]->name() << "\"" << std::endl;
            }
            continue;
        }

        // Skip observation if we have a binned observation
        if (obs->eventtype() == "CountsCube") {
            if (logTerse()) {
                log << "Warning: Skipping binned observation \"";
                log << obs->name()+"\"" << std::endl;
            }
            continue;
        }

        // Write header for current observation
        if (logTerse()) {
            if (obs->name().length() > 1) {
                log.header3("Observation "+obs->name());
            }
            else {
                log.header3("Observation");
            }
        }

        // Map events into sky map
        map_events(obs);
            
        // Increment observation counter
        num_obs++;

    } // endfor: looped over observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save observation
 *
 * This method saves the counts map(s) into (a) FITS file(s).
 ***************************************************************************/
void ctskymap::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save sky map");
    }

    // Get output filename
    m_outmap = (*this)["outmap"].filename();

    // Save sky map
    m_skymap.save(m_outmap, clobber());

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
void ctskymap::get_parameters(void)
{
    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {

        // Throw exception if no input observation file is given
        require_inobs(G_GET_PARAMETERS);

        // Get observation container without response (not needed)
        m_obs = get_observations(false);

    } // endif: there was no observation in the container

    // Create sky map based on task parameters
    m_skymap = create_map(m_obs);

    // Get remaining parameters
    m_emin = (*this)["emin"].real();
    m_emax = (*this)["emax"].real();

    // Read ahead parameters
    if (read_ahead()) {
        m_outmap = (*this)["outmap"].filename();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Map events into a sky map
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::no_list
 *            No event list found in observation.
 * @exception GCTAException::no_pointing
 *            No valid CTA pointing found.
 *
 * This method maps the events found in a CTA events list into a sky map.
 ***************************************************************************/
void ctskymap::map_events(GCTAObservation* obs)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Make sure that the observation holds a CTA event list. If this
        // is not the case then throw an exception.
        if (dynamic_cast<const GCTAEventList*>(obs->events()) == NULL) {
            throw GException::no_list(G_BIN_EVENTS);
        }

        // Setup energy range covered by data
        GEnergy  emin;
        GEnergy  emax;
        GEbounds ebds;
        emin.TeV(m_emin);
        emax.TeV(m_emax);

        // Initialise binning statistics
        int num_outside_map    = 0;
        int num_outside_erange = 0;
        int num_in_map         = 0;

        // Fill sky map
        GCTAEventList* events = static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));
        for (int i = 0; i < events->size(); ++i) {

            // Get event
            GCTAEventAtom* event = (*events)[i];

            // Skip if energy is out of range
            if (event->energy() < emin || event->energy() > emax) {
                num_outside_erange++;
                continue;
            }

            // Determine sky pixel
            GCTAInstDir* inst  = (GCTAInstDir*)&(event->dir());
            GSkyDir      dir   = inst->dir();
            GSkyPixel    pixel = m_skymap.dir2pix(dir);

            // Skip if pixel is out of range
            if (pixel.x() < -0.5 || pixel.x() > (m_skymap.nx() - 0.5) ||
                pixel.y() < -0.5 || pixel.y() > (m_skymap.ny() - 0.5)) {
                num_outside_map++;
                continue;
            }

            // Fill event in skymap
            m_skymap(pixel, 0) += 1.0;
            num_in_map++;

        } // endfor: looped over all events

        // Log binning results
        if (logTerse()) {
            log << std::endl;
            log.header1("Mapping");
            log << gammalib::parformat("Events in list");
            log << obs->events()->size() << std::endl;
            log << gammalib::parformat("Events in map");
            log << num_in_map << std::endl;
            log << gammalib::parformat("Events outside map area");
            log << num_outside_map << std::endl;
            log << gammalib::parformat("Events outside energy range");
            log << num_outside_erange << std::endl;
        }

        // Log map
        if (logTerse()) {
            log << std::endl;
            log.header1("Sky map");
            log << m_skymap << std::endl;
        }

    } // endif: observation was valid

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
void ctskymap::init_members(void)
{
    // Initialise members
    m_obs.clear();
    m_skymap.clear();
    m_emin = 0.0;
    m_emax = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctskymap::copy_members(const ctskymap& app)
{
    // Copy attributes
    m_obs    = app.m_obs;
    m_skymap = app.m_skymap;
    m_emin   = app.m_emin;
    m_emax   = app.m_emax;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctskymap::free_members(void)
{
    // Return
    return;
}
