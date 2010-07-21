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
 * @brief Setup Crab point source powerlaw model.
 ***************************************************************************/
GModels crab_plaw(void)
{
    // Setup GModels for optimizing
    GModelSpatialPtsrc point_source;
    GModelSpectralPlaw power_law;
    GModel             crab;
    GModels            models;
    try {
        GSkyDir dir;
        //dir.radec_deg(83.6331, +22.0145);
        dir.radec_deg(117.0, -33.0);  // Adapt to source position in file
        point_source = GModelSpatialPtsrc(dir);
        power_law    = GModelSpectralPlaw(1.0e-7, -2.1);
        power_law.par(0)->min(1.0e-12);
        crab         = GModel(point_source, power_law);
        crab.name("Crab");
        models.append(crab);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable setup GModels for optimizing." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }

    // Return model
    return models;
}


/***********************************************************************//**
 * @brief Run maximum likelihood analysis
 ***************************************************************************/
int ctlike::run(void)
{
    // Initialise return code
    int rc = 0;
    
    // Branch on method
    std::string method = par("method")->value();
    if (toupper(method) == "UNBINNED")
        rc = unbinned();
    else
        rc = binned();
    
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
    std::string caldb  = par("caldb")->value();
    std::string irf    = par("irf")->value();

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
    run.response(irf,caldb);
    run.roi(roi);
    run.gti()->add(tstart, tstop);
    run.ebounds()->append(emin, emax);
    obs.append(run);

    // Setup GModels for optimizing
    GModels models = crab_plaw();
    obs.models(models);

    // Perform LM optimization
    GOptimizerLM opt;
    opt.max_iter(1000);
    obs.optimize(opt);
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
    std::string caldb  = par("caldb")->value();
    std::string irf    = par("irf")->value();

    // Declare observations
    GObservations   obs;
    GCTAObservation run;

    // Load binned CTA observation
    run.load_binned(cntmap);
    run.response(irf,caldb);
    obs.append(run);

    // Setup GModels for optimizing
    GModels models = crab_plaw();
    obs.models(models);

    // Perform LM optimization
    GOptimizerLM opt;
    opt.max_iter(1000);
    obs.optimize(opt);
    std::cout << obs << std::endl;

    // Return
    return rc;
}
