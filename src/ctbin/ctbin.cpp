/***************************************************************************
 *                      ctbin - CTA data binning tool                      *
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


/***********************************************************************//**
 * @brief Main entry point of application
 *
 * @param[in] argc Number of command line arguments.
 * @param[in] argv Command line arguments.
 *
 * This is the main entry point of the gtbin application. It allocates a
 * ctbin object and runs the application.
 ***************************************************************************/
int main (int argc, char *argv[])
{
    // Create instance of application
    ctbin application(argc, argv);
    
    // Run application
    int rc = application.run();

    // Return
    return rc;
}


/***********************************************************************//**
 * @brief Run gtbin application
 ***************************************************************************/
int ctbin::run(void)
{
    // Test dump
    std::cout << *this << std::endl;

    // Get events and output file names
    std::string evfile  = par("evfile")->value();
    std::string outfile = par("outfile")->value();
    
    // Get clobber
    std::string sclobber = par("clobber")->value();
    int clobber = 0;
    if (tolower(sclobber) == "yes")
        clobber = 1;

    // Load events as unbinned data
    GCTAObservation obs;
    obs.load_unbinned(evfile);

    // Setup energy range covered by data
    GEnergy emin;
    GEnergy emax;
    emin.TeV(todouble(par("emin")->value()));
    emax.TeV(todouble(par("emax")->value()));
    int enumbins = toint(par("enumbins")->value());
    obs.ebounds()->setlog(emin, emax, enumbins);

    // Setup time range covered by data
    GTime tstart;
    GTime tstop;
    tstart.met(0.0);
    tstop.met(10000.0);
    obs.gti()->add(tstart, tstop);

    std::cout << obs;
    
    // Create skymap
    std::string wcs    = par("proj")->value();
    std::string coords = par("coordsys")->value();
    double      x      = todouble(par("xref")->value());
    double      y      = todouble(par("yref")->value());
    double      dx     = todouble(par("binsz")->value());
    double      dy     = todouble(par("binsz")->value());
    int         nx     = toint(par("nxpix")->value());
    int         ny     = toint(par("nypix")->value());
    GSkymap map = GSkymap(wcs, coords, x, y, dx, dy, nx, ny, enumbins);

    // Fill sky map
    GCTAEventList* events = (GCTAEventList*)obs.events();
    for (GCTAEventList::iterator event = events->begin(); event != events->end(); ++event) {

        // Determine sky pixel
        GSkyDir   dir   = ((GCTAInstDir*)event->dir())->skydir();
        GSkyPixel pixel = map.dir2xy(dir);

        // Skip if pixel is out of range
        if (pixel.x() < -0.5 || pixel.x() > (nx-0.5) ||
            pixel.y() < -0.5 || pixel.y() > (ny-0.5))
            continue;

        // Determine energy bin. Skip if we are outside the energy range
        int index = obs.ebounds()->index(*(event->energy()));
        if (index == -1)
            continue;

        // Fill event in skymap
        map(pixel, index) += 1.0;
        
    } // endfor: looped over all events

    std::cout << map;
    
    // Save counts map
    GFits file;
    map.write(&file);
    obs.ebounds()->write(&file);
    obs.gti()->write(&file);
    file.saveto(outfile, clobber);

    // Return
    return 0;
}