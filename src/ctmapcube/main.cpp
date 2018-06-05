/***************************************************************************
 *                 ctmapcube - CTA map cube generation tool                *
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
 * @file ctmapcube/main.cpp
 * @brief CTA map cube generation tool main code
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "support.hpp"
#include "ctmapcube.hpp"


/***********************************************************************//**
 * @brief Main entry point of application
 *
 * @param[in] argc Number of command line arguments.
 * @param[in] argv Command line arguments.
 *
 * This is the main entry point of the ctmapcube application. It allocates a
 * ctmapcube object and executes the application. Any exceptions that occur
 * will be catched and corresponding error messages written in the
 * application logger and into the standard output.
 ***************************************************************************/
int main (int argc, char *argv[])
{
    // Initialise return code
    int rc = 1;

    // Initialise pointer on application
    ctmapcube* application = NULL;

    // Execute application
    try {

        // Create instance of application
        application = new ctmapcube(argc, argv);

        // Execute ctool
        rc = execute_ctool(application);

    }

    catch (std::exception &e) {

        // Report exception
        report_ctool_failure("ctmapcube", e.what());

    }

    // Delete application
    delete application;

    // Return
    return rc;
}
