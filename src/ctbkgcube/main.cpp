/***************************************************************************
 *                  ctbkgcube - CTA background cube tool                   *
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
 * @file ctbkgcube/main.cpp
 * @brief CTA background cube tool main code
 * @author Chia-Chun Lu
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "support.hpp"
#include "ctbkgcube.hpp"


/***********************************************************************//**
 * @brief Main entry point of application
 *
 * @param[in] argc Number of command line arguments.
 * @param[in] argv Command line arguments.
 *
 * This is the main entry point of the ctbkgcube application. It allocates a
 * ctbkgcube object and executes the application. Any exceptions that occur
 * will be catched and corresponding error messages written in the
 * application logger and into the standard output.
 ***************************************************************************/
int main (int argc, char *argv[])
{
    // Initialise return code
    int rc = 1;

    // Initialise pointer on application
    ctbkgcube* application = NULL;

    // Execute application
    try {

        // Create instance of application
        application = new ctbkgcube(argc, argv);

        // Execute ctool
        rc = execute_ctool(application);

    }

    catch (std::exception &e) {

        // Report exception
        report_ctool_failure("ctbkgcube", e.what());

    }

    // Delete application
    delete application;

    // Return
    return rc;
}
