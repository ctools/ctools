/***************************************************************************
 *                 ctexpcube - Exposure cube generation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Juergen Knoedlseder                         *
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

    // Set energy dispersion flag for all CTA observations and save old
    // values in save_edisp vector
    std::vector<bool> save_edisp;
    save_edisp.assign(m_obs.size(), false);
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            save_edisp[i] = obs->response()->apply_edisp();
            obs->response()->apply_edisp(m_apply_edisp);
        }
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
        log.header1("Generate exposure cube");
    }

    // Fill exposure
    m_expcube.fill(m_obs);

    // Restore energy dispersion flag for all CTA observations
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            obs->response()->apply_edisp(save_edisp[i]);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save exposure cube
 ***************************************************************************/
void ctexpcube::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save exposure cube");
    }

    // Get output filename
    m_outcube = (*this)["outcube"].filename();

    // Save exposure cube
    m_expcube.save(m_outcube, clobber());

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
    // Initialise members
    m_outcube.clear();
    m_apply_edisp = false;

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
    // Copy attributes
    m_outcube     = app.m_outcube;
    m_apply_edisp = app.m_apply_edisp;

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

        // Build observation container
        m_obs = get_observations();

    } // endif: there was no observation in the container

    // Read energy dispersion flag
    m_apply_edisp = (*this)["edisp"].boolean();

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

    // Read output filename (if needed)
    if (read_ahead()) {
        m_outcube = (*this)["outcube"].filename();
    }

    // Return
    return;
}
