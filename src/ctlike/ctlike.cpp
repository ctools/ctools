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
 * @file ctlike.cpp
 * @brief CTA maximum likelihood tool implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include "ctlike.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Run maximum likelihood analysis
 ***************************************************************************/
int ctlike::run(void)
{
    // Initialise return code
    int rc = 0;

    // Get standard parameters
    m_method = toupper(par("method")->value());
    m_caldb  = par("caldb")->value();
    m_irf    = par("irf")->value();
    m_srcmdl = par("srcmdl")->value();
    m_models = GModels(m_srcmdl);

    // Set fixed parameters
    m_max_iter = 1000;
    m_opt      = new GOptimizerLM;

    // Set optimizer parameters
    ((GOptimizerLM*)m_opt)->max_iter(m_max_iter);

    // Branch on method
    if (m_method == "UNBINNED")
        rc = unbinned();
    else
        rc = binned();

    // Free optimizer
    delete m_opt;

    // Return
    return rc;
}


/***********************************************************************//**
 * @brief Perform unbinned maximum likelihood analysis
 ***************************************************************************/
int ctlike::unbinned(void)
{
    // Initialise return code
    int rc = 0;

    // Get parameters
    std::string evfile = par("evfile")->value();

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Setup ROI covered by data
    GCTAInstDir instDir;
    GCTARoi     roi;
    instDir.radec_deg(117.02, -33.35);  // Adapt to file
    roi.centre(instDir);
    roi.radius(2.5);

    // Setup energy range covered by data
    GEnergy emin;
    GEnergy emax;
    emin.TeV(0.02);
    emax.TeV(100.0);

    // Setup time range covered by data
    GTime tstart;
    GTime tstop;
    tstart.met(0.0);
    tstop.met(10000.0);

    // Load data and response and set ROI, energy range and time range
    // for analysis
    run.load_unbinned(evfile);
    run.response(m_irf, m_caldb);
    run.roi(roi);
    run.gti()->add(tstart, tstop);
    run.ebounds()->append(emin, emax);
    obs.append(run);

    // Setup models for optimizing
    obs.models(m_models);

    // Perform LM optimization
    obs.optimize(*m_opt);
    std::cout << obs << std::endl;
    
    // Return
    return rc;
}


/***********************************************************************//**
 * @brief Perform binned maximum likelihood analysis
 ***************************************************************************/
int ctlike::binned(void)
{
    // Initialise return code
    int rc = 0;

    // Get parameters
    std::string cntmap = par("cntmap")->value();

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load binned CTA observation
    run.load_binned(cntmap);
    run.response(m_irf, m_caldb);
    obs.append(run);

    // Setup models for optimizing
    obs.models(m_models);

    // Perform LM optimization
    obs.optimize(*m_opt);
    std::cout << obs << std::endl;

    // Return
    return rc;
}
