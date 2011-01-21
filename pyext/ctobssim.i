/***************************************************************************
 *               ctobssim - CTA observation simulation tool                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file ctobssim.i
 * @brief CTA observation simulation tool python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctobssim.hpp"
%}
%include gammalib.i


/***********************************************************************//**
 * @class ctobssim
 *
 * @brief CTA data selection tool python interface defintion.
 ***************************************************************************/
class ctobssim : public GApplication  {
public:
    // Constructors and destructors
    ctobssim(void);
    ctobssim(int argc, char *argv[]);
    ~ctobssim(void);

    // Methods
    void           run(void);
    void           get_parameters(void);
    GPhotons       simulate_photons(void);
    GCTAEventList* simulate_events(const GPhotons& photons);
    void           save_events(const GCTAEventList* events);
};


/***********************************************************************//**
 * @brief CTA observation simulation tool extension
 ***************************************************************************/
%extend ctobssim {
    ctobssim copy() {
        return (*self);
    }
}
