/***************************************************************************
 *                        ctool_obs - [WHAT] tool                          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file ctool_obs.cpp
 * @brief [WHAT] tool implementation
 * @author [AUTHOR]
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ctool_obs.hpp"
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __ OpenMP section _____________________________________________________ */
//#ifdef _OPENMP
//#include <omp.h>
//#endif

/* __ Method name definitions ____________________________________________ */

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
 * Constructs empty [what] tool.
 ***************************************************************************/
ctool_obs::ctool_obs(void) : ctobservation(CTOOL_OBS_NAME, VERSION)
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
 * Constructs [what] tool from an observation container.
 ***************************************************************************/
ctool_obs::ctool_obs(const GObservations& obs) : ctobservation(CTOOL_OBS_NAME, VERSION, obs)
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
 * Constructs [what] tool using command line arguments for user
 * parameter setting.
 ***************************************************************************/
ctool_obs::ctool_obs(int argc, char *argv[]) : ctobservation(CTOOL_OBS_NAME, VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app [WHAT] tool.
 *
 * Constructs [what] tool from another [what] tool.
 ***************************************************************************/
ctool_obs::ctool_obs(const ctool_obs& app) : ctobservation(app)
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
 * Destructs [what] tool.
 ***************************************************************************/
ctool_obs::~ctool_obs(void)
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
 * @param[in] app [WHAT] tool.
 * @return [WHAT] tool.
 *
 * Assigns [what] tool.
 ***************************************************************************/
ctool_obs& ctool_obs::operator=(const ctool_obs& app)
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
 * @brief Clear [WHAT] tool
 *
 * Clears [what] tool.
 ***************************************************************************/
void ctool_obs::clear(void)
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
 * @brief Process [what] tool
 ***************************************************************************/
void ctool_obs::process(void)
{
    // Get task parameters
    get_parameters();

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // TODO: Your code goes here

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save something
 *
 * Saves something.
 ***************************************************************************/
void ctool_obs::save(void)
{
    // Write header
    log_header1(TERSE, "Save something");

    // TODO: Your code goes here

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
void ctool_obs::init_members(void)
{
    // Initialise members
    // TODO: Your code goes here

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app [WHAT] tool.
 ***************************************************************************/
void ctool_obs::copy_members(const ctool_obs& app)
{
    // Copy attributes
    // TODO: Your code goes here

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctool_obs::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * @todo Implement method
 ***************************************************************************/
void ctool_obs::get_parameters(void)
{
    // TODO: Your code goes here

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}
