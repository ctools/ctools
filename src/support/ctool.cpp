/***************************************************************************
 *                        ctool - ctool base class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file ctool.hpp
 * @brief ctool base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ctool.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_EBOUNDS                                  "ctool::get_ebounds()"

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
ctool::ctool(void) : GApplication()
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Name constructor
 *
 * @param[in] name Application name.
 * @param[in] version Application version.
 ***************************************************************************/
ctool::ctool(const std::string& name, const std::string& version) :
       GApplication(name, version)
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Command line constructor
 *
 * @param[in] name Application name.
 * @param[in] version Application version.
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 ***************************************************************************/
ctool::ctool(const std::string& name, const std::string& version,
             int argc, char *argv[]) : 
       GApplication(name, version, argc, argv)
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Application.
 ***************************************************************************/
ctool::ctool(const ctool& app) : GApplication(app)
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
ctool::~ctool(void)
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
ctool& ctool::operator=(const ctool& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->GApplication::operator=(app);

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
void ctool::init_members(void)
{
    // Initialise members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctool::copy_members(const ctool& app)
{
    // Copy attributes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctool::free_members(void)
{
    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get the energy boundaries
 *
 * @exception GException::invalid_value
 *            No valid energy boundary extension found.
 *
 * Get the energy boundaries according to the user parameters. The method
 * supports loading of energy boundary information from the EBOUNDS or
 * ENERGYBINS extension, or setting energy boundaries using a linear
 * or logarithmical spacing.
 ***************************************************************************/
GEbounds ctool::get_ebounds(void)
{
    // Allocate energy boundaries
    GEbounds ebounds;

    // Get energy binning algorithm
    std::string ebinalg = (*this)["ebinalg"].string();

    // If energy binning algorithm is of type "FILE" (case sensitive), then
    // read energy boundaries from FITS file ...
    if (ebinalg == "FILE") {

        // Get filename
        std::string ebinfile = (*this)["ebinfile"].filename();

        // Open energy boundary file using the EBOUNDS or ENERGYBINS
        // extension. Throw an exception if opening fails.
        GFits file(ebinfile);
        if (file.contains("EBOUNDS")) {
            file.close();
            ebounds.load(ebinfile,"EBOUNDS");
        }
        else if (file.contains("ENERGYBINS")) {
            file.close();
            ebounds.load(ebinfile,"ENERGYBINS");
        }
        else {
            file.close();
            std::string msg = "No extension with name \"EBOUNDS\" or"
                              " \"ENERGYBINS\" found in FITS file"
                              " \""+ebinfile+"\".\n"
                              "An \"EBOUNDS\" or \"ENERGYBINS\" extension"
                              " is required if the parameter \"ebinalg\""
                              " is set to \"FILE\".";
            throw GException::invalid_value(G_GET_EBOUNDS, msg);
        }
        
    } // endif: ebinalg was "FILE"

    // ... otherwise use a linear or a logarithmically-spaced energy binning
    else {

        // Get task parameters
    	double emin     = (*this)["emin"].real();
    	double emax     = (*this)["emax"].real();
    	int    enumbins = (*this)["enumbins"].integer();

        // Initialise log mode for ebinning
        bool log = true;

        // Check if algorithm is linear
        if (ebinalg == "LIN") {
            log = false;
        }

        // Setup energy bins
        ebounds = GEbounds(enumbins, GEnergy(emin, "TeV"),
                                     GEnergy(emax, "TeV"), log);

    } // endelse: ebinalg was not "FILE"

    // Return energy boundaries
    return ebounds;
}
