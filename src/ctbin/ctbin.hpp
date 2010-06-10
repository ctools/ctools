/***************************************************************************
 *                      ctbin - CTA data binning tool                      *
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
 * @file ctbin.hpp
 * @brief CTA data binning tool definition
 * @author J. Knodlseder
 */

#ifndef CTBIN_HPP
#define CTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define APP_NAME    "ctbin"
#define APP_VERSION "v1r0p0"


/***********************************************************************//**
 * @class ctbin
 *
 * @brief CTA data binning tool interface defintion.
 ***************************************************************************/
class ctbin : public GApplication  {
public:
    // Constructors and destructors
    ctbin(int argc, char *argv[]) : 
                  GApplication(APP_NAME, APP_VERSION, argc, argv) { return; }

    // Methods
    int run(void);
};

#endif /* CTBIN_HPP */
