/***************************************************************************
 *                 ctexpcube - Exposure cube generation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Juergen Knoedlseder                         *
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
 * @file ctexpcube.cpp
 * @brief Exposure cube generation tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctexpcube.hpp"
#include "GTools.hpp"
#include "GWcs.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                        "ctexpcube::get_parameters()"
#define G_SET_FROM_CNTMAP          "ctexpcube::set_from_cntmap(std::string&)"

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
ctexpcube::ctexpcube(void) : ctool(CTEXPCUBE_NAME, CTEXPCUBE_VERSION)
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
ctexpcube::ctexpcube(const GObservations& obs) :
           ctool(CTEXPCUBE_NAME, CTEXPCUBE_VERSION)
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
ctexpcube::ctexpcube(int argc, char *argv[]) :
           ctool(CTEXPCUBE_NAME, CTEXPCUBE_VERSION, argc, argv)
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
ctexpcube::ctexpcube(const ctexpcube& app) : ctool(app)
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
ctexpcube::~ctexpcube(void)
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
ctexpcube& ctexpcube::operator=(const ctexpcube& app)
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
void ctexpcube::clear(void)
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
 * @brief Generate the exposure cube(s).
 *
 * This method reads the task parameters from the parfile, sets up the
 * observation container, loops over all CTA observations in the container
 * and generates an exposure cube from the CTA observations.
 ***************************************************************************/
void ctexpcube::run(void)
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

    // Warn if there are not enough energy bins
    int    n          = m_expcube.energies().size();
    double logEmin    = std::log10(m_expcube.energies()[0].TeV());
    double logEmax    = std::log10(m_expcube.energies()[n-1].TeV());
    int    nrequired  = int((logEmax - logEmin) * 25.0);
    if ((n < nrequired) && logTerse()) {
        log << std::endl;
        log << "WARNING: ";
        log << "Only " << n-1 << " energy bins have been requested. This may be "
               "too few energy bins.";
        log << std::endl << "         ";
        log << "At least 25 bins per decade in energy are recommended, which "
               "in your case";
        log << std::endl << "         ";
        log << "would be " << nrequired-1 << " bins. Consider increasing the "
               "number of energy bins.";
        log << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Initialise exposure cube");
    }

    // Initialise exposure cube
    init_cube();

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        log.header1(gammalib::number("Observation",m_obs.size()));
        log << m_obs.print(m_chatter) << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Generate exposure cube");
    }

    // Fill exposure cube
    m_expcube.fill(m_obs, &log);

    // Log exposure cube
    if (logTerse()) {
        log << m_expcube.print(m_chatter) << std::endl;
    }

    // Optionally publish exposure cube
    if (m_publish) {
        publish();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save exposure cube
 *
 * Saves the exposure cube into a FITS file.
 ***************************************************************************/
void ctexpcube::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save exposure cube");
    }

    // Get exposure cube filename
    m_outcube = (*this)["outcube"].filename();

    // Save only if filename is non-empty
    if (!m_outcube.is_empty()) {

        // Log filename
        if (logTerse()) {
            log << "Save exposure cube into file \""+m_outcube+"\".";
            log << std::endl;
        }

        // Save exposure cube
        m_expcube.save(m_outcube, clobber());

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Publish exposure cube
 *
 * @param[in] name Exposure cube name.
 ***************************************************************************/
void ctexpcube::publish(const std::string& name)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Publish exposure cube");
    }

    // Set default name is user name is empty
    std::string user_name(name);
    if (user_name.empty()) {
        user_name = CTEXPCUBE_NAME;
    }

    // Log filename
    if (logTerse()) {
        log << "Publish \""+user_name+"\" exposure cube." << std::endl;
    }

    // Publish exposure cube
    m_expcube.cube().publish(user_name);

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
void ctexpcube::init_members(void)
{
    // Initialise user parameters
    m_outcube.clear();
    m_addbounds = true;
    m_publish   = false;
    m_chatter   = static_cast<GChatter>(2);

    // Initialise protected members
    m_obs.clear();
    m_expcube.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctexpcube::copy_members(const ctexpcube& app)
{
    // Copy user parameters
    m_outcube   = app.m_outcube;
    m_addbounds = app.m_addbounds;
    m_publish   = app.m_publish;
    m_chatter   = app.m_chatter;

    // Copy protected members
    m_obs     = app.m_obs;
    m_expcube = app.m_expcube;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctexpcube::free_members(void)
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
void ctexpcube::get_parameters(void)
{
    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {

        // Throw exception if no input observation file is given
        require_inobs(G_GET_PARAMETERS);

        // Throw exception if counts cube is given
        require_inobs_nocube(G_GET_PARAMETERS);

        // Build observation container
        m_obs = get_observations();

    } // endif: there was no observation in the container

    // Get the incube filename
    std::string incube = (*this)["incube"].filename();

    // Check for filename validity
    if ((incube == "NONE") || (gammalib::strip_whitespace(incube) == "")) {

       // Create an event cube based on task parameters
       GCTAEventCube cube = create_cube(m_obs);

       // Define exposure cube
       m_expcube = GCTACubeExposure(cube);

    } // endif: filename was not valid

    // ... otherwise setup the exposure cube from the counts map
    else {

        // Load event cube from filename
        GCTAEventCube cube(incube);

        // Define exposure cube
        m_expcube = GCTACubeExposure(cube);

    } // endelse: cube loaded from file

    // Get remaining parameters
    m_addbounds = (*this)["addbounds"].boolean();
    m_publish   = (*this)["publish"].boolean();
    m_chatter   = static_cast<GChatter>((*this)["chatter"].integer());

    // Read output filename (if needed)
    if (read_ahead()) {
        m_outcube = (*this)["outcube"].filename();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise exposure cube
 *
 * Initialise the exposure cube.
 ***************************************************************************/
void ctexpcube::init_cube(void)
{
    // Extract exposure cube definition
    const GWcs* proj   = static_cast<const GWcs*>(m_expcube.cube().projection());
    std::string wcs    = m_expcube.cube().projection()->code();
    std::string coords = m_expcube.cube().projection()->coordsys();
    double      x      = proj->crval(0);
    double      y      = proj->crval(1);
    double      dx     = proj->cdelt(0);
    double      dy     = proj->cdelt(1);
    int         nx     = m_expcube.cube().nx();
    int         ny     = m_expcube.cube().ny();

    // Extract energies
    GEnergies energies = m_expcube.energies();

    // If requested, insert energies at all event list energy boundaries
    if (m_addbounds) {

        // Set logger
        GLog* logger = NULL;
        if (logTerse()) {
            logger = &log;
        }
    
        // Loop over all observations
        for (int i = 0; i < m_obs.size(); ++i) {
    
            // Get observation and continue only if it is a CTA observation
            const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>
                                         (m_obs[i]);

            // Skip observation if it's not a CTA observation
            if (cta == NULL) {
                continue;
            }

            // Skip observation if it does not contain an event list
            if (cta->eventtype() != "EventList") {
                continue;
            }

            // Insert energy boundaries
            energies = insert_energy_boundaries(energies, *cta, logger);

        } // endfor: looped over all observations

       } // endif: energy bin insertion requested

    // Setup exposure cube
    m_expcube = GCTACubeExposure(wcs, coords, x, y, dx, dy, nx, ny, energies);

    // Log exposure cube
    if (logTerse()) {
        log << m_expcube.print(m_chatter) << std::endl;
    }

    // Return
    return;
};
