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
 * @brief CTA observation simulation tool Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctobssim.hpp"
#include "GTools.hpp"
%}
%include gammalib.i


/***********************************************************************//**
 * @class ctobssim
 *
 * @brief CTA data selection tool Python interface
 ***************************************************************************/
class ctobssim : public GApplication  {
public:
    // Constructors and destructors
    ctobssim(void);
    explicit ctobssim(GObservations obs);
    ctobssim(int argc, char *argv[]);
    ctobssim(const ctobssim& app);
    virtual ~ctobssim(void);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void);
    void           get_parameters(void);
    void           set_list(GCTAObservation* obs);
    GPhotons       simulate_photons(const GCTAObservation* obs, const GModels& models);
    void           simulate_source(GCTAObservation* obs, const GPhotons& photons);
    void           simulate_background(GCTAObservation* obs, const GModels& models);
};


/***********************************************************************//**
 * @brief CTA observation simulation tool Python extension
 ***************************************************************************/
%extend ctobssim {
    ctobssim copy() {
        return (*self);
    }
}
