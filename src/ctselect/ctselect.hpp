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
 * @file ctselect.hpp
 * @brief CTA data selection tool definition
 * @author J. Knodlseder
 */

#ifndef CTSELECT_HPP
#define CTSELECT_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"

/* __Definitions _________________________________________________________ */
#define APP_NAME    "ctselect"
#define APP_VERSION "v1r0p0"


/***********************************************************************//**
 * @class ctselect
 *
 * @brief CTA data selection tool interface defintion.
 ***************************************************************************/
class ctselect : public GApplication  {
public:
    // Constructors and destructors
    ctselect(int argc, char *argv[]) : 
                  GApplication(APP_NAME, APP_VERSION, argc, argv) { return; }

    // Methods
    int run(void);
};

#endif /* CTSELECT_HPP */
