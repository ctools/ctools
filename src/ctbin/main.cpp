/***************************************************************************
 *                    ctbin - CTA data binning main code                   *
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
 * @file main.cpp
 * @brief CTA data binning tool main code
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include "ctbin.hpp"


/***********************************************************************//**
 * @brief Main entry point of application
 *
 * @param[in] argc Number of command line arguments.
 * @param[in] argv Command line arguments.
 *
 * This is the main entry point of the gtbin application. It allocates a
 * ctbin object and runs the application.
 ***************************************************************************/
int main (int argc, char *argv[])
{
    // Create instance of application
    ctbin application(argc, argv);
    
    // Run application
    int rc = application.run();

    // Return
    return rc;
}
