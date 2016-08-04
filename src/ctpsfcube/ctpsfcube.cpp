/***************************************************************************
 *                  ctpsfcube - PSF cube generation tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Chia-Chun Lu                                *
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
 * @file ctpsfcube.cpp
 * @brief PSF cube generation tool implementation
 * @author Chia-Chun Lu
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctpsfcube.hpp"
#include "GTools.hpp"
#include "GWcs.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                        "ctpsfcube::get_parameters()"

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
ctpsfcube::ctpsfcube(void) : ctool(CTPSFCUBE_NAME, CTPSFCUBE_VERSION)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observations constructor
 *
 * @param[in] obs Observation container.
 *
 * This method creates an instance of the class by copying an existing
 * observations container.
 ***************************************************************************/
ctpsfcube::ctpsfcube(const GObservations& obs) :
           ctool(CTPSFCUBE_NAME, CTPSFCUBE_VERSION)
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
ctpsfcube::ctpsfcube(int argc, char *argv[]) :
           ctool(CTPSFCUBE_NAME, CTPSFCUBE_VERSION, argc, argv)
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
ctpsfcube::ctpsfcube(const ctpsfcube& app) : ctool(app)
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
ctpsfcube::~ctpsfcube(void)
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
ctpsfcube& ctpsfcube::operator=(const ctpsfcube& app)
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
 * @brief Clear ctpsfcube tool
 *
 * Clears ctpsfcube tool.
 ***************************************************************************/
void ctpsfcube::clear(void)
{
    // Free members
    free_members();
    this->ctool::free_members();

    // Clear base class (needed to conserve tool name and version)
    this->GApplication::clear();

    // Initialise members
    this->ctool::init_members();
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
 * and generates a PSF cube from the CTA observations.
 ***************************************************************************/
void ctpsfcube::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Initialise PSF cube");
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
        log.header1("Generate PSF cube");
    }

    // Fill PSF cube
    m_psfcube.fill(m_obs, &log);

    // Log PSF cube
    if (logTerse()) {
        log << m_psfcube.print(m_chatter) << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save PSF cube
 ***************************************************************************/
void ctpsfcube::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save PSF cube");
    }

    // Get output filename
    m_outcube = (*this)["outcube"].filename();

    // Save only if filename is non-empty
    if (!m_outcube.is_empty()) {

        // Log filename
        if (logTerse()) {
            log << gammalib::parformat("PSF cube file");
            log << m_outcube.url() << std::endl;
        }

        // Save PSF cube
        m_psfcube.save(m_outcube, clobber());

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
void ctpsfcube::init_members(void)
{
    // Initialise user parameters
    m_outcube.clear();
    m_addbounds = true;
    m_chatter   = static_cast<GChatter>(2);

    // Initialise protected members
    m_obs.clear();
    m_psfcube.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctpsfcube::copy_members(const ctpsfcube& app)
{
    // Copy user parameters
    m_outcube   = app.m_outcube;
    m_addbounds = app.m_addbounds;
    m_chatter   = app.m_chatter;

    // Copy protected members
    m_obs     = app.m_obs;
    m_psfcube = app.m_psfcube;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctpsfcube::free_members(void)
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
void ctpsfcube::get_parameters(void)
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

    // Get additional binning parameters
    double amax     = (*this)["amax"].real();
    int    anumbins = (*this)["anumbins"].integer();

    // Check for filename validity
    if ((incube == "NONE") || (gammalib::strip_whitespace(incube) == "")) {

        // Create an event cube based on task parameters
        GCTAEventCube cube = create_cube(m_obs);

        // Define psf cube
        m_psfcube = GCTACubePsf(cube, amax, anumbins);

    }

    // ... otherwise setup the exposure cube from the counts map
    else {

        // Load event cube from filename
        GCTAEventCube cube(incube);

        // Define psf cube
        m_psfcube = GCTACubePsf(cube, amax, anumbins);

    } // endelse: cube loaded from file

    // Get remaining parameters
    m_addbounds = (*this)["addbounds"].boolean();
    m_chatter   = static_cast<GChatter>((*this)["chatter"].integer());

    // Read output filename (if needed)
    if (read_ahead()) {
        m_outcube = (*this)["outcube"].filename();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise PSF cube
 *
 * Initialise the PSF cube.
 ***************************************************************************/
void ctpsfcube::init_cube(void)
{
    // Extract exposure cube definition
    const GWcs* proj   = static_cast<const GWcs*>(m_psfcube.cube().projection());
    std::string wcs    = m_psfcube.cube().projection()->code();
    std::string coords = m_psfcube.cube().projection()->coordsys();
    double      x      = proj->crval(0);
    double      y      = proj->crval(1);
    double      dx     = proj->cdelt(0);
    double      dy     = proj->cdelt(1);
    int         nx     = m_psfcube.cube().nx();
    int         ny     = m_psfcube.cube().ny();
    double      dmax   = (*this)["amax"].real();
    int         ndbins = (*this)["anumbins"].integer();

    // Extract energies
    GEnergies energies = m_psfcube.energies();

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

    // Setup PSF cube
    m_psfcube = GCTACubePsf(wcs, coords, x, y, dx, dy, nx, ny, energies,
                            dmax, ndbins);

    // Log PSF cube
    if (logTerse()) {
        log << m_psfcube.print(m_chatter) << std::endl;
    }

    // Return
    return;
};
