/***************************************************************************
 *           ctobssim - CTA observation simulator tool main code           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file ctobssim.cpp
 * @brief CTA observation simulator tool main code
 * @author J. Knodlseder
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
 ***************************************************************************/
int main (int argc, char *argv[])
{
    // Initialise return code
    int rc = 1;

    // Create instance of application
    ctobssim application(argc, argv);

    // Run application
    try {
        // Run application
        application.run();

        // Save results
        application.save();

        // Signal success
        rc = 0;
    }
    catch (std::exception &e) {

        // Extract error message
        std::string message = e.what();
        std::string signal  = "*** ERROR encounterted in the execution of"
                              " ctobssim. Run aborted ...";

        // Write error in logger
        application.log << signal  << std::endl;
        application.log << message << std::endl;

        // Write error on standard output
        std::cout << signal  << std::endl;
        std::cout << message << std::endl;

    } // endcatch: catched any application error

    // Return
    return rc;
}
