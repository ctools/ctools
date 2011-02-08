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
#define CTOBSSIM_VERSION "00-02-00"


/***********************************************************************//**
 * @class ctobssim
 *
 * @brief CTA observation simulator tool interface defintion
 *
 * This class simulates CTA observation(s) using Monte Carlo sampling of the
 * source and background models. The class supports simulation of data of
 * multiple CTA observations in one shot. If multiple CTA observations are
 * processed and the save method is called, events FITS files will be written
 * for each observation.
 ***************************************************************************/
class ctobssim : public GApplication  {
public:
    // Constructors and destructors
    ctobssim(void);
    explicit ctobssim(GObservations obs);
    ctobssim(int argc, char *argv[]);
    ctobssim(const ctobssim& app);
    virtual ~ctobssim(void);

    // Operators
    ctobssim& operator= (const ctobssim& app);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           set_list(GCTAObservation* obs);
    GPhotons       simulate_photons(const GCTAObservation* obs, const GModels& models);
    void           simulate_source(GCTAObservation* obs, const GPhotons& photons);
    void           simulate_background(GCTAObservation* obs, const GModels& models);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctobssim& app);
    void free_members(void);

    // User parameters
    std::string   m_infile;     //!< Input model
    std::string   m_outfile;    //!< Output events file
    std::string   m_caldb;      //!< Calibration database repository
    std::string   m_irf;        //!< Instrument response function
    int           m_seed;       //!< Random number generator seed 
    double        m_ra;         //!< RA of pointing direction
    double        m_dec;        //!< DEC of pointing direction
    double        m_rad;        //!< FOV radius
    double        m_tmin;       //!< Start time (MET)
    double        m_tmax;       //!< Stop time (MET)
    double        m_emin;       //!< Lower energy (TeV)
    double        m_emax;       //!< Upper energy (TeV)

    // Protected members
    double        m_area;       //!< Surface area for simulation (cm2)
    GRan          m_ran;        //!< Random number generator
    GObservations m_obs;        //!< Observation container
};

#endif /* CTOBSSIM_HPP */
