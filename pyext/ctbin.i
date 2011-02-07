/***************************************************************************
 *                      ctbin - CTA data binning tool                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file ctbin.i
 * @brief CTA data binning tool Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctbin.hpp"
%}
%include gammalib.i


/***********************************************************************//**
 * @class ctbin
 *
 * @brief CTA data binning tool Python interface defintion
 ***************************************************************************/
class ctbin : public GApplication  {
public:
    // Constructors and destructors
    ctbin(void);
    explicit ctbin(GObservations obs);
    ctbin(int argc, char *argv[]);
    virtual ~ctbin(void);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           bin_events(GCTAObservation* obs);
};


/***********************************************************************//**
 * @brief CTA data binning tool SWIG extension
 ***************************************************************************/
%extend ctbin {
    ctbin copy() {
        return (*self);
    }
}
