/***************************************************************************
 *                    ctselect - CTA data selection tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file ctselect.cpp
 * @brief CTA data selection tool implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include "ctselect.hpp"


/***********************************************************************//**
 * @brief Main entry point
 *
 * @param[in] argc Number of arguments
 * @param[in] argv Arguments
 ***************************************************************************/
int main (int argc, char *argv[])
{
    // Initialise return code
    int rc = 0;

    // Put all code in try-catch brackets to catch all errors
    try {

        // Create instance of application
        ctselect application(argc, argv);

        // Run application
        rc = application.run();

    }
    catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        rc = -1;
    }

    // Return
    return rc;
}


/***********************************************************************//**
 * @brief Application
 ***************************************************************************/
int ctselect::run(void)
{
    // Test dump
    std::cout << *this << std::endl;

    // DUMMY FOR TESTING
    std::string value;
    value = par("string")->value();
    std::cout << value << std::endl;
    value = par("chatter")->value();
    std::cout << value << std::endl;
    value = par("clobber")->value();
    std::cout << value << std::endl;

    // Return
    return 0;
}
