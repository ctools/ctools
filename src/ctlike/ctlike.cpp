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
 * @file ctlike.cpp
 * @brief CTA maximum likelihood tool implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
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
ctlike::ctlike(void) : GApplication(CTLIKE_NAME, CTLIKE_VERSION)
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

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

    // Write header into logger
    log_header();

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
 *
 * The following analysis steps are performed:
 * 1. Read the parameters (and write them into logger)
 * 2. Load observation
 * 3. Setup models for optimizing
 * 4. Optimize model (and write result into logger)
 * 5. Write results out as XML model
 ***************************************************************************/
void ctlike::run(void)
{
    // Switch screen logging on in debug mode
    if (logDebug())
        log.cout(true);

    // Get parameters
    get_parameters();

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Load observation
    if (m_method == "UNBINNED")
        unbinned(m_evfile);
    else
        binned(m_cntmap);

    // Write observation into logger
    if (logTerse()) {
        log << std::endl;
        log.header1("Observation");
        log << m_obs << std::endl;
    }

    // Setup models for optimizing.
    m_models = GModels(m_srcmdl);
    m_obs.models(m_models);

    // Optimize
    optimize_lm();

    // Write results into logger
    if (logTerse()) {
        log << " Maximum log likelihood ....: " << m_logL << std::endl;
        log << " Npred .....................: " << m_obs.npred() << std::endl;
        log << m_models << std::endl << std::endl;
    }
    
    // Write results out as XML model
    if (toupper(m_outmdl) != "NONE")
        m_models.save(m_outmdl);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user.
 ***************************************************************************/
void ctlike::get_parameters(void)
{
    // Get standard parameters
    m_method = toupper(par("method")->value());
    m_stat   = toupper(par("stat")->value());
    m_refit  = par("refit")->boolean();
    m_caldb  = par("caldb")->value();
    m_irf    = par("irf")->value();
    m_srcmdl = par("srcmdl")->value();
    m_outmdl = par("outmdl")->value();

    // Get unbinned parameters ...
    if (m_method == "UNBINNED") {
        m_evfile = par("evfile")->value();
    }
    
    // ... or get binned parameters
    else if (m_method == "BINNED") {
        m_cntmap = par("cntmap")->value();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Optimize using Levenberg-Marquardt method
 ***************************************************************************/
void ctlike::optimize_lm(void)
{
    // Allocate optimizer. The logger is only passed to the optimizer
    // constructor if optimizer logging is requested.
    GOptimizerLM* opt = (logTerse()) ? new GOptimizerLM(log) 
                                     : new GOptimizerLM();
    
    // Assign optimizer
    m_opt = opt;

    // Set optimizer parameters
    opt->max_iter(m_max_iter);

    // Write Header for optimization and indent for optimizer logging
    if (logTerse()) {
        log << std::endl;
        log.header1("Maximum likelihood optimisation");
        log.indent(1);
    }
    
    // Perform LM optimization
    m_obs.optimize(*m_opt);

    // Optionally refit
    if (m_refit)
        m_obs.optimize(*m_opt);

    // Get models back
    m_models = *(m_obs.models());

    // Store maximum log likelihood value
    m_logL = -opt->value();

    // Write optimization results
    log.indent(0);
    if (logTerse()) {
        log << std::endl;
        log.header1("Maximum likelihood optimisation results");
        log << *opt << std::endl;
    }

    // Free optimizer
    delete opt;
    
    // Reset optimizer pointer
    m_opt = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load unbinned observation and append to observations contained
 *
 * @param[in] evfile Events file name.
 ***************************************************************************/
void ctlike::unbinned(const std::string& evfile)
{
    // Declare CTA observation
    GCTAObservation obs;

    // Set statistics
    obs.statistics(m_stat);

    // Load data
    obs.load_unbinned(evfile);

    // Set reponse
    obs.response(m_irf, m_caldb);

    // Append observation to contained
    m_obs.append(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load binned observation and append to observations contained
 *
 * @param[in] cntmap Counts map file name.
 ***************************************************************************/
void ctlike::binned(const std::string& cntmap)
{
    // Declare CTA observation
    GCTAObservation obs;

    // Set statistics
    obs.statistics(m_stat);

    // Load data
    obs.load_binned(cntmap);

    // Set reponse
    obs.response(m_irf, m_caldb);

    // Append observation to contained
    m_obs.append(obs);

    // Return
    return;
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
    m_obs.clear();
    m_max_iter = 100;   // Set maximum number of iterations to 100
    m_logL     = 0.0;
    m_opt      = NULL;

    // Set logger properties
    log.date(true);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctlike::free_members(void)
{
    // Free members

    // Return
    return;
}
