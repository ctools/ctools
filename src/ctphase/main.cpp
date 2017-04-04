/***************************************************************************
 *                ctphase - Append phase information to CTA events file    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Leonardo Di Venere                         *
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
 * @file ctphase/main.cpp
 * @brief Append phase information to CTA events file
 * @author Leonardo Di Venere
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "support.hpp"
#include "ctphase.hpp"


/***********************************************************************//**
 * @brief Main entry point
 *
 * @param[in] argc Number of arguments
 * @param[in] argv Arguments
 *
 * This is the main entry point of the ctphase application. It allocates a
 * ctphase object and executes the application. Any exceptions that occur
 * will be catched and corresponding error messages written in the
 * application logger and into the standard output.
 ***************************************************************************/
int main (int argc, char *argv[])
{
    // Initialise return code
    int rc = 1;

    // Initialise pointer on application
    ctphase* application = NULL;

    // Execute application
    try {

        // Create instance of application
        application = new ctphase(argc, argv);

        // Execute ctool
        rc = execute_ctool(application);

    }

    catch (std::exception &e) {

        // Report exception
        report_ctool_failure("ctphase", e.what());

    }

    // Delete application
    delete application;

    // Return
    return rc;
}
