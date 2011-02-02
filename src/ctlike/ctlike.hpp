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
 * @file ctlike.hpp
 * @brief CTA maximum likelihood tool definition
 * @author J. Knodlseder
 */

#ifndef CTLIKE_HPP
#define CTLIKE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTLIKE_NAME    "ctlike"
#define CTLIKE_VERSION "00-01-00"


/***********************************************************************//**
 * @class ctlike
 *
 * @brief CTA maximum likelihood tool interface defintion.
 ***************************************************************************/
class ctlike : public GApplication  {
public:
    // Constructors and destructors
    ctlike(void);
    ctlike(int argc, char *argv[]);
    ~ctlike(void);

    // Methods
    void run(void);
    void get_parameters(void);
    void optimize_lm(void);
    void unbinned(const std::string& evfile);
    void binned(const std::string& cntmap);

protected:
    // Protected methods
    void init_members(void);
    void free_members(void);

    // User parameters
    std::string   m_method;     //!< Likelihood fitting method (binned/unbinned)
    std::string   m_stat;       //!< Optimisation statistics (poisson/gaussian)
    bool          m_refit;      //!< Refitting
    std::string   m_caldb;      //!< Calibration database
    std::string   m_irf;        //!< Instrument response functions
    std::string   m_srcmdl;     //!< Source model XML file
    std::string   m_outmdl;     //!< Source model output XML file
    GModels       m_models;     //!< Source models
    GObservations m_obs;        //!< Observations

    // Unbinned likelihood parameters
    std::string   m_evfile;     //!< Events file

    // Binned likelihood parameters
    std::string   m_cntmap;     //!< Counts map

    // Internal parameters
    int           m_max_iter;   //!< Maximum number of iterations
    double        m_logL;       //!< Maximum log likelihood
    GOptimizer*   m_opt;        //!< Optimizer
};

#endif /* CTLIKE_HPP */
