/***************************************************************************
 *                   ctlike - CTA maximum likelihood tool                  *
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
 * @file ctlike.i
 * @brief CTA maximum likelihood tool Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctlike.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class ctlike
 *
 * @brief CTA maximum likelihood tool Python interface
 ***************************************************************************/
class ctlike : public GApplication  {
public:
    // Constructors and destructors
    ctlike(void);
    explicit ctlike(GObservations obs);
    ctlike(int argc, char *argv[]);
    ctlike(const ctlike& app);
    virtual ~ctlike(void);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           optimize_lm(void);
    void           load_unbinned(const std::string& evfile);
    void           load_binned(const std::string& cntmap);
};


/***********************************************************************//**
 * @brief CTA maximum likelihood tool Python extension
 ***************************************************************************/
%extend ctlike {
    ctlike copy() {
        return (*self);
    }
}
