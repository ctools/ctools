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
#define CTLIKE_VERSION "v1r0p0"


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
    int run(void);
    int unbinned(void);
    int binned(void);

protected:
    // Protected methods
    void init_members(void);
    void free_members(void);

    // User parameters
    std::string   m_method;     //!< Likelihood fitting method (binned/unbinned)
    std::string   m_caldb;      //!< Calibration database
    std::string   m_irf;        //!< Instrument response functions
    std::string   m_srcmdl;     //!< Source model XML file
    std::string   m_outmdl;     //!< Source model output XML file
    GModels       m_models;     //!< Source models
    GObservations m_obs;        //!< Observations

    // Fixed parameters
    int           m_max_iter;   //!< Maximum number of iterations
    GOptimizer*   m_opt;        //!< Optimizer
};

#endif /* CTLIKE_HPP */
