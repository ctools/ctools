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
 * @file ctselect.i
 * @brief CTA data selection tool SWIG definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctselect.hpp"
%}
%include gammalib.i


/***********************************************************************//**
 * @class ctselect
 *
 * @brief CTA data selection tool SWIG interface defintion.
 ***************************************************************************/
class ctselect : public GApplication  {
public:
    // Constructors and destructors
    ctselect(void);
    ctselect(int argc, char *argv[]);
    ~ctselect(void);

    // Methods
    void run(void);
    void get_parameters(void);
    void select(void);
    void append_gti(void);
    void copy(void);
};


/***********************************************************************//**
 * @brief CTA data selection tool SWIG extension
 ***************************************************************************/
%extend ctselect {
    ctselect copy() {
        return (*self);
    }
}
