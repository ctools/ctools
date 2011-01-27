/***************************************************************************
 *                ctobssim - CTA observation simulator tool                *
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
 * @file ctobssim.hpp
 * @brief CTA observation simulator tool definition
 * @author J. Knodlseder
 */

#ifndef CTOBSSIM_HPP
#define CTOBSSIM_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTOBSSIM_NAME    "ctobssim"
#define CTOBSSIM_VERSION "v1r0p0"


/***********************************************************************//**
 * @class ctobssim
 *
 * @brief CTA observation simulator tool interface defintion.
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
    void           set_events_keywords(GFits* file);

protected:
    // Protected methods
    void init_members(void);
    void free_members(void);

    // User parameters
    std::string     m_infile;     //!< Input model
    std::string     m_outfile;    //!< Output events file
    std::string     m_caldb;      //!< Calibration database repository
    std::string     m_irf;        //!< Instrument response function
    int             m_seed;       //!< Random number generator seed 
    double          m_ra;         //!< RA of pointing direction
    double          m_dec;        //!< DEC of pointing direction
    double          m_rad;        //!< FOV radius
    double          m_tmin;       //!< Start time (MET)
    double          m_tmax;       //!< Stop time (MET)
    double          m_emin;       //!< Lower energy (TeV)
    double          m_emax;       //!< Upper energy (TeV)
    double          m_area;       //!< Surface area for simulation (cm2)
    GRan            m_ran;        //!< Random number generator
    GCTAPointing    m_pnt;        //!< CTA pointing
    GCTAObservation m_obs;        //!< CTA observation
    GModels         m_models;     //!< Models
};

#endif /* CTOBSSIM_HPP */
