/***************************************************************************
 *                      ctbin - CTA data binning tool                      *
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
 * @file ctbin.cpp
 * @brief CTA data binning tool implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include "ctbin.hpp"
#include "GTools.hpp"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
ctbin::ctbin(void) : GApplication(CTBIN_NAME, CTBIN_VERSION)
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
ctbin::ctbin(int argc, char *argv[]) : 
                          GApplication(CTBIN_NAME, CTBIN_VERSION, argc, argv)
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
ctbin::~ctbin(void)
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
 * @brief Run gtbin application
 ***************************************************************************/
void ctbin::run(void)
{
    // Get parameters
    get_parameters();

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Bin data
    bin();

    // Write separator into logger
    if (logTerse())
        log << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user.
 ***************************************************************************/
void ctbin::get_parameters(void)
{
    // Get parameters
    m_evfile   = par("evfile")->value();
    m_outfile  = par("outfile")->value();
    m_emin     = par("emin")->real();
    m_emax     = par("emax")->real();
    m_enumbins = par("enumbins")->integer();
    m_proj     = par("proj")->value();
    m_coordsys = par("coordsys")->value();
    m_xref     = par("xref")->real();
    m_yref     = par("yref")->real();
    m_binsz    = par("binsz")->real();
    m_nxpix    = par("nxpix")->integer();
    m_nypix    = par("nypix")->integer();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Bin the data
 ***************************************************************************/
void ctbin::bin(void)
{
    // Load events as unbinned data
    GCTAObservation obs;
    obs.load_unbinned(m_evfile);

    // Setup energy range covered by data
    GEnergy  emin;
    GEnergy  emax;
    GEbounds ebds;
    emin.TeV(m_emin);
    emax.TeV(m_emax);
    ebds.setlog(emin, emax, m_enumbins);
    obs.ebounds(ebds);

    // Log observation
    if (logTerse()) {
        log << std::endl;
        log.header1("Observation");
        log << obs << std::endl;
    }

    // Create skymap
    GSkymap map = GSkymap(m_proj, m_coordsys,
                          m_xref, m_yref, m_binsz, m_binsz,
                          m_nxpix, m_nypix, m_enumbins);

    // Initialise binning statistics
    int num_outside_map  = 0;
    int num_outside_ebds = 0;
    int num_in_map       = 0;
    
    // Fill sky map
    GCTAEventList* events = (GCTAEventList*)obs.events();
    for (GCTAEventList::iterator event = events->begin(); event != events->end(); ++event) {

        // Determine sky pixel
        GCTAInstDir* inst  = (GCTAInstDir*)&(event->dir());
        GSkyDir      dir   = inst->skydir();
        GSkyPixel    pixel = map.dir2xy(dir);

        // Skip if pixel is out of range
        if (pixel.x() < -0.5 || pixel.x() > (m_nxpix-0.5) ||
            pixel.y() < -0.5 || pixel.y() > (m_nypix-0.5)) {
            num_outside_map++;
            continue;
        }

        // Determine energy bin. Skip if we are outside the energy range
        int index = obs.ebounds().index(event->energy());
        if (index == -1) {
            num_outside_ebds++;
            continue;
        }

        // Fill event in skymap
        map(pixel, index) += 1.0;
        num_in_map++;
        
    } // endfor: looped over all events

    // Log binning results
    if (logTerse()) {
        log << std::endl;
        log.header1("Binning");
        log << parformat("Events in list");
        log << obs.events()->size() << std::endl;
        log << parformat("Events in map");
        log << num_in_map << std::endl;
        log << parformat("Events outside map area");
        log << num_outside_map << std::endl;
        log << parformat("Events outside energy bins");
        log << num_outside_ebds << std::endl;
    }

    // Log observation
    if (logTerse()) {
        log << std::endl;
        log.header1("Counts map");
        log << map << std::endl;
    }
    
    // Save counts map
    GFits file;
    map.write(&file);
    obs.ebounds().write(&file);
    obs.gti().write(&file);
    file.saveto(m_outfile, clobber());

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
void ctbin::init_members(void)
{
    // Initialise members
    m_evfile.clear();
    m_outfile.clear();
    m_emin     = 0.0;
    m_emax     = 0.0;
    m_enumbins = 0;
    m_proj.clear();
    m_coordsys.clear();
    m_xref  = 0.0;
    m_yref  = 0.0;
    m_binsz = 0.0;
    m_nxpix = 0;
    m_nypix = 0;

    // Set logger properties
    log.date(true);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctbin::free_members(void)
{
    // Free members

    // Return
    return;
}
