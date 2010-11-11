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


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
ctlike::ctlike(void) :  GApplication(CTLIKE_NAME, CTLIKE_VERSION)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Command line constructor
 *
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 ***************************************************************************/
ctlike::ctlike(int argc, char *argv[]) : 
                        GApplication(CTLIKE_NAME, CTLIKE_VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
ctlike::~ctlike(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

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
    m_outmdl = par("outmdl")->value();
    m_models = GModels(m_srcmdl);

    // Set fixed parameters
    m_opt    = new GOptimizerLM;

    // Set optimizer parameters
    ((GOptimizerLM*)m_opt)->max_iter(m_max_iter);

    // Branch on method
    if (m_method == "UNBINNED")
        rc = unbinned();
    else
        rc = binned();

    // Write result as XML model
    if (toupper(m_outmdl) != "NONE")
        m_models.save(m_outmdl);

    // Dump results
    std::cout << m_models << std::endl;

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
    m_obs.append(run);

    // Setup models for optimizing
    m_obs.models(m_models);

    // Perform LM optimization
    m_obs.optimize(*m_opt);
    
    // Get models back
    m_models = *(m_obs.models());
    
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
    GCTAObservation run;

    // Load binned CTA observation
    run.load_binned(cntmap);
    run.response(m_irf, m_caldb);
    m_obs.append(run);

    // Setup models for optimizing
    m_obs.models(m_models);

    // Perform LM optimization
    m_obs.optimize(*m_opt);

    // Get models back
    m_models = *(m_obs.models());

    // Return
    return rc;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void ctlike::init_members(void)
{
    // Initialise members
    m_method.clear();
    m_caldb.clear();
    m_irf.clear();
    m_srcmdl.clear();
    m_outmdl.clear();
    m_models.clear();
    //m_obs.clear();     //!< NOT YET IMPLEMENTED
    m_max_iter = 1000;
    m_opt      = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctlike::free_members(void)
{
    // Free members
    delete m_opt;

    // Return
    return;
}
