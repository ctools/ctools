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
ctedispcube::ctedispcube(void) : ctobservation(CTEDISPCUBE_NAME, CTEDISPCUBE_VERSION)
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
             ctobservation(CTEDISPCUBE_NAME, CTEDISPCUBE_VERSION, obs)
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
 * Constructs an energy dispersion tool by parsing the arguments provided
 * on the command line.
 ***************************************************************************/
ctedispcube::ctedispcube(int argc, char *argv[]) :
             ctobservation(CTEDISPCUBE_NAME, CTEDISPCUBE_VERSION, argc, argv)
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
ctedispcube::ctedispcube(const ctedispcube& app) : ctobservation(app)
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
 * @brief Clear energy dispersion tool
 *
 * Set the energy disperison tool to an empty tool.
 ***************************************************************************/
void ctedispcube::clear(void)
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

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Initialise energy dispersion cube
    init_cube();

    // Write header into logger
    log_header1(TERSE, "Generate energy dispersion cube");

    // Set pointer to logger dependent on chattiness
    GLog* logger = (logNormal()) ? &log : NULL;

    // Fill energy dispersion cube
    m_edispcube.fill(m_obs, logger);

    // Write energy dispersion cube into logger
    log_string(EXPLICIT, m_edispcube.print(m_chatter));

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
    // Write header into logger
    log_header1(TERSE, "Save energy dispersion cube");

    // Get output filename
    m_outcube = (*this)["outcube"].filename();

    // Save energy dispersion cube if filename and the energy dispersion cube
    // are not empty
    if (!m_outcube.is_empty() && !m_edispcube.cube().is_empty()) {
        m_edispcube.save(m_outcube, clobber());
    }

    // Write into logger what has been done
    std::string fname = (m_outcube.is_empty()) ? "NONE" : m_outcube.url();
    if (m_edispcube.cube().is_empty()) {
        fname.append(" (cube is empty, no file created)");
    }
    log_value(NORMAL, "Energy dispersion cube file", fname);

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
    // Initialise user parameters
    m_outcube.clear();
    m_addbounds = true;
    m_chatter   = static_cast<GChatter>(2);

    // Initialise protected members
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
    // Copy user parameters
    m_outcube   = app.m_outcube;
    m_addbounds = app.m_addbounds;
    m_chatter   = app.m_chatter;

    // Copy protected members
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
    // Setup observations from "inobs" parameter. Do not accept counts cubes.
    setup_observations(m_obs, true, true, false);

    // Get the incube filename
    std::string incube = (*this)["incube"].filename();

    // Get additional binning parameters
    double migramax  = (*this)["migramax"].real();
    int    migrabins = (*this)["migrabins"].integer();

    // If the "incube" file name is valid then setup the energy dispersion
    // cube from the counts cube. Otherwise create a counts cube from the
    // user parameters
    GCTAEventCube cube = is_valid_filename(incube) ? GCTAEventCube(incube)
                                                   : create_cube(m_obs);

    // Define edisp cube
    m_edispcube = GCTACubeEdisp(cube, migramax, migrabins);

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
 * @brief Initialise energy dispersion cube
 *
 * Initialise the energy dispersion cube.
 ***************************************************************************/
void ctedispcube::init_cube(void)
{
    // Write header into logger
    log_header1(TERSE, "Initialise energy dispersion cube");

    // Extract energy dispersion definition
    const GWcs* proj   = static_cast<const GWcs*>(m_edispcube.cube().projection());
    std::string wcs    = m_edispcube.cube().projection()->code();
    std::string coords = m_edispcube.cube().projection()->coordsys();
    double      x      = proj->crval(0);
    double      y      = proj->crval(1);
    double      dx     = proj->cdelt(0);
    double      dy     = proj->cdelt(1);
    int         nx     = m_edispcube.cube().nx();
    int         ny     = m_edispcube.cube().ny();
    double      mmax   = (*this)["migramax"].real();
    int         nmbins = (*this)["migrabins"].integer();

    // Extract energies
    GEnergies energies = m_edispcube.energies();

    // If requested, insert energies at all event list energy boundaries
    if (m_addbounds) {

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
            energies = insert_energy_boundaries(energies, *cta);

        } // endfor: looped over all observations

       } // endif: energy bin insertion requested

    // Setup energy dispersion cube
    m_edispcube = GCTACubeEdisp(wcs, coords, x, y, dx, dy, nx, ny, energies,
                                mmax, nmbins);

    // Write energy dispersion cube in logger
    log_string(NORMAL, m_edispcube.print(m_chatter));

    // Return
    return;
}
