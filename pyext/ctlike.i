/***************************************************************************
 *                   ctlike - CTA maximum likelihood tool                  *
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
 * @file ctlike.i
 * @brief CTA maximum likelihood tool SWIG definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctlike.hpp"
%}
%include gammalib.i


/***********************************************************************//**
 * @class ctlike
 *
 * @brief CTA maximum likelihood tool SWIG interface defintion.
 ***************************************************************************/
class ctlike : public GApplication  {
public:
    // Constructors and destructors
    ctlike(void) : GApplication(CTLIKE_NAME, CTLIKE_VERSION) { return; }

    // Methods
    int run(void);
};


/***********************************************************************//**
 * @brief CTA maximum likelihood tool SWIG extension
 ***************************************************************************/
%extend ctlike {
    ctlike copy() {
        return (*self);
    }
}
