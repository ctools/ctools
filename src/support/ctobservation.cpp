/***************************************************************************
 *             ctobservation - Base class for observation tools            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file ctobservation.cpp
 * @brief Observation tool base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ctobservation.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Name constructor
 *
 * @param[in] name Observation tool name.
 * @param[in] version Observation tool version.
 *
 * Constructs a observation tool from the @p name and @p version. See the
 * equivalent ctool constructor for details.
 ***************************************************************************/
ctobservation::ctobservation(const std::string& name,
                             const std::string& version) : ctool(name, version)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observations constructor
 *
 * @param[in] name Observation tool name.
 * @param[in] version Observation tool version.
 * param[in] obs Observation container.
 *
 * Constructs a observation tool from the @p name, @p version and an
 * observation container.
 ***************************************************************************/
ctobservation::ctobservation(const std::string&   name,
                             const std::string&   version,
                             const GObservations& obs) : ctool(name, version)
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
 * @param[in] name Observation tool name.
 * @param[in] version Observation tool version.
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 *
 * Constructs a observation tool from the @p name, @p version and command
 * line arguments. See the equivalent ctool constructor for details.
 ***************************************************************************/
ctobservation::ctobservation(const std::string& name,
                             const std::string& version,
                             int   argc,
                             char *argv[]) : ctool(name, version, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Observation tool.
 *
 * Constructs an instance of a observation tool by copying information from
 * another observation tool.
 ***************************************************************************/
ctobservation::ctobservation(const ctobservation& app) : ctool(app)
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
 * Destructs the observation tool.
 ***************************************************************************/
ctobservation::~ctobservation(void)
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
 * @param[in] app Observation tool.
 * @return Observation tool.
 *
 * Assigns a observation tool.
 ***************************************************************************/
ctobservation& ctobservation::operator=(const ctobservation& app)
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


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void ctobservation::init_members(void)
{
    // Initialise protected members
    m_obs.clear();

    // Initialise private members
    m_index_unbinned = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Observation tool.
 ***************************************************************************/
void ctobservation::copy_members(const ctobservation& app)
{
    // Copy protected members
    m_obs = app.m_obs;

    // Copy private members
    m_index_unbinned = app.m_index_unbinned;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctobservation::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return first unbinned CTA observation (const version)
 *
 * @return Const pointer to first unbinned CTA observation
 *
 * Returns a const pointer to the first unbinned CTA observation in the
 * container. If no CTA observation exists a NULL pointer is returned.
 *
 * The method calls next_unbinned_observation(). See the method for details.
 ***************************************************************************/
const GCTAObservation* ctobservation::first_unbinned_observation(void) const
{
    // Initialise index
    m_index_unbinned = 0;

    // Get next unbinned CTA observation
    const GCTAObservation* obs = next_unbinned_observation();

    // Return first CTA observation
    return obs;
}


/***********************************************************************//**
 * @brief Return next unbinned CTA observation (const version)
 *
 * @return Const pointer to next unbinned CTA observation
 *
 * Returns a const pointer to the next unbinned CTA observation in the
 * container. If no CTA observation exists any more a NULL pointer is
 * returned.
 *
 * The method writes for each encountered observation a level 3 header into
 * the logger. It will also signal when an observation was skipped because
 * it either was not a CTA observation or not an unbinned observation.
 *
 * @todo Logger methods should be declared const to avoid the const casting.
 ***************************************************************************/
const GCTAObservation* ctobservation::next_unbinned_observation(void) const
{
    // Initialise pointer on CTA observation
    const GCTAObservation* obs = NULL;

    // Loop over all remaining observation in the container
    for (; m_index_unbinned < m_obs.size(); ++m_index_unbinned) {

        // Write header for the current observation
        const_cast<ctobservation*>(this)->log_header3(TERSE,
                                   get_obs_header(m_obs[m_index_unbinned]));

        // Case the observation to a CTA observation. This will return a
        // NULL pointer if the observation is not a CTA observation.
        obs = dynamic_cast<const GCTAObservation*>(m_obs[m_index_unbinned]);

        // Skip observation if it's not CTA
        if (obs == NULL) {
            std::string msg = " Skipping "+
                              m_obs[m_index_unbinned]->instrument()+
                              " observation";
            const_cast<ctobservation*>(this)->log_string(NORMAL, msg);
            continue;
        }

        // Skip observation if we have a binned observation
        if (obs->eventtype() == "CountsCube") {
            obs             = NULL;
            std::string msg = " Skipping binned "+
                              m_obs[m_index_unbinned]->instrument()+
                              " observation";
            const_cast<ctobservation*>(this)->log_string(NORMAL, msg);
            continue;
        }

        // If we come to this point we have an unbinned CTA observation and
        // we can forward the index to the next index and break the loop
        m_index_unbinned++;
        break;

    } // endfor: looped over all observation

    // Return next CTA observation
    return obs;
}
