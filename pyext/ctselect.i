/***************************************************************************
 *                    ctselect - CTA data selection tool                   *
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
 * @file ctselect.i
 * @brief CTA data selection tool Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctselect.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class ctselect
 *
 * @brief CTA data selection tool Python interface
 ***************************************************************************/
class ctselect : public GApplication  {
public:
    // Constructors and destructors
    ctselect(void);
    explicit ctselect(GObservations obs);
    ctselect(int argc, char *argv[]);
    ctselect(const ctselect& app);
    virtual ~ctselect(void);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           select_events(GCTAObservation* obs, const std::string& filename);
};


/***********************************************************************//**
 * @brief CTA data selection tool Python extension
 ***************************************************************************/
%extend ctselect {
    ctselect copy() {
        return (*self);
    }
}
