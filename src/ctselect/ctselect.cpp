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
    // Create instance of application
    ctselect application(argc, argv);
    
    // Run application
    int rc = application.run();

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
    
    // Return
    return 0;
}