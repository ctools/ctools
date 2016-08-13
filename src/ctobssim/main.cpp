/***************************************************************************
 *           ctobssim - CTA observation simulator tool main code           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2016 by Juergen Knoedlseder                         *
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
 * @file ctobssim/main.cpp
 * @brief CTA observation simulator tool main code
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctobssim.hpp"


/***********************************************************************//**
 * @brief Main entry point
 *
 * @param[in] argc Number of arguments.
 * @param[in] argv Arguments.
 *
 * This is the main entry point of the ctobssim application. It allocates a
 * ctobssim object and executes the application. Any exceptions that occur
 * will be catched and corresponding error messages written in the
 * application logger and on the standard output.
 ***************************************************************************/
int main (int argc, char *argv[])
{
    // Initialise return code
    int rc = 1;

    // Initialise pointer on application
    ctobssim* application = NULL;

    // Execute application
    try {

        // Create instance of application
        application = new ctobssim(argc, argv);

        // Execute application
        application->execute();

        // Delete application
        delete application;

        // Signal success
        rc = 0;
    }
    catch (std::exception &e) {

        // Extract error message
        std::string message = e.what();
        std::string signal  = "*** ERROR encounterted in the execution of "
                              "ctobssim. Run aborted ...";

        // If application exists then write error in logger and delete
        // application
        if (application != NULL) {
            application->log << signal  << std::endl;
            application->log << message << std::endl;
            delete application;
        }

        // Write error on standard output
        std::cout << signal  << std::endl;
        std::cout << message << std::endl;

    } // endcatch: catched any application error

    // Return
    return rc;
}
