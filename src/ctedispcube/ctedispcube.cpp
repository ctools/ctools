/***************************************************************************
 *          ctedispcube - Energy dispersion cube generation tool           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Maria Haupt                                      *
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
 * @file ctedispcube.cpp
 * @brief Energy dispersion cube generation tool implementation
 * @author Maria Haupt
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctedispcube.hpp"
#include "GTools.hpp"
#include "GWcs.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                      "ctedispcube::get_parameters()"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty energy dispersion tool.
 ***************************************************************************/
ctedispcube::ctedispcube(void) : ctool(CTEDISPCUBE_NAME, CTEDISPCUBE_VERSION)
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
 * Constructs an energy dispersion tool from the information that is provided
 * in an observation container @p obs.
 ***************************************************************************/
ctedispcube::ctedispcube(const GObservations& obs) :
             ctool(CTEDISPCUBE_NAME, CTEDISPCUBE_VERSION)
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
 * Constructs an energy dispersion tool by parsing the arguments provided
 * on the command line.
 ***************************************************************************/
ctedispcube::ctedispcube(int argc, char *argv[]) :
             ctool(CTEDISPCUBE_NAME, CTEDISPCUBE_VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Energy dispersion tool.
 *
 * Constructs an energy dispersion tool by copying anothere energy
 * dispersion tool.
 ***************************************************************************/
ctedispcube::ctedispcube(const ctedispcube& app) : ctool(app)
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
 * Desctructs an energy dispersion tool.
 ***************************************************************************/
ctedispcube::~ctedispcube(void)
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
 * @param[in] app Energy dispersion tool.
 * @return Energy dispersion tool.
 *
 * Assigns energy dispersion tool.
 ***************************************************************************/
ctedispcube& ctedispcube::operator=(const ctedispcube& app)
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
 * @brief Clear energy dispersion tool
 *
 * Set the energy disperison tool to an empty tool.
 ***************************************************************************/
void ctedispcube::clear(void)
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
 * @brief Generate the energy dispersion cube
 *
 * Generates the energy dispersion cube by looping over all unbinned CTA
 * observations in the observation container.
 ***************************************************************************/
void ctedispcube::run(void)
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
        log.header1("Generate energy dispersion cube");
    }

    // Fill Edisp
    m_edispcube.fill(m_obs, &log);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save energy dispersion cube
 *
 * Saves the energy dispersion cube in a FITS file. The FITS filename is
 * provided by the "outcube" parameter.
 ***************************************************************************/
void ctedispcube::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save Edisp cube");
    }

    // Get output filename
    m_outcube = (*this)["outcube"].filename();

    // Save only if filename is non-empty
    if (!m_outcube.is_empty()) {

        // Log filename
        if (logTerse()) {
            log << "Save energy dispersion cube into file \""+m_outcube+"\".";
            log << std::endl;
        }

        // Save energy dispersion cube
        m_edispcube.save(m_outcube, clobber());

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
void ctedispcube::init_members(void)
{
    // Initialise members
    m_outcube.clear();
    m_obs.clear();
    m_edispcube.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Energy dispersion cube.
 ***************************************************************************/
void ctedispcube::copy_members(const ctedispcube& app)
{
    // Copy members
    m_outcube   = app.m_outcube;
    m_obs       = app.m_obs;
    m_edispcube = app.m_edispcube;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctedispcube::free_members(void)
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
void ctedispcube::get_parameters(void)
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
    double migramax  = (*this)["migramax"].real();
    int    migrabins = (*this)["migrabins"].integer();

    // Check for filename validity
    if ((incube == "NONE") || (gammalib::strip_whitespace(incube) == "")) {

        // Create an event cube based on task parameters
        GCTAEventCube cube = create_cube(m_obs);

        // Define edisp cube
        m_edispcube = GCTACubeEdisp(cube, migramax, migrabins);

    }

    // ... otherwise setup the exposure cube from the counts map
    else {

        // Load event cube from filename
        GCTAEventCube cube(incube);

        // Define Edisp cube
        m_edispcube = GCTACubeEdisp(cube, migramax, migrabins);

    } // endelse: cube loaded from file

    // Read output filename (if needed)
    if (read_ahead()) {
        m_outcube = (*this)["outcube"].filename();
    }

    // Return
    return;
}
