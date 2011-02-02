/***************************************************************************
 *                    ctselect - CTA data selection tool                   *
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
 * @file ctselect.cpp
 * @brief CTA data selection tool implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include "ctselect.hpp"
#include "GTools.hpp"

/* __ Debug definitions __________________________________________________ */
#define BIT_COLUMN_TEST    0
#define STRING_COLUMN_TEST 1

/* __ Coding definitions _________________________________________________ */
#define SELECT_TIME   1
#define SELECT_ENERGY 1
#define SELECT_ROI    1


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
ctselect::ctselect(void) : GApplication(CTSELECT_NAME, CTSELECT_VERSION)
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
ctselect::ctselect(int argc, char *argv[]) : 
                    GApplication(CTSELECT_NAME, CTSELECT_VERSION, argc, argv)
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
ctselect::~ctselect(void)
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
 * @brief Select event data
 *
 * A GTI is appended to the FITS file as in the actual data format no GTI
 * exists.
 ***************************************************************************/
void ctselect::run(void)
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

    // Select events
    select();

    // Append GTI to FITS file
    append_gti();

    // Set data selection keywords
    write_ds_keys();

    // Save output FITS file
    m_file.saveto(m_outfile, clobber());

    // Close FITS file
    m_file.close();

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
void ctselect::get_parameters(void)
{
    // Get parameters
    m_infile  = par("infile")->value();
    m_outfile = par("outfile")->value();
    m_ra      = par("ra")->real();
    m_dec     = par("dec")->real();
    m_rad     = par("rad")->real();
    m_tmin    = par("tmin")->real();
    m_tmax    = par("tmax")->real();
    m_emin    = par("emin")->real();
    m_emax    = par("emax")->real();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Select events
 *
 * Select events from a FITS file by making use of the selection possibility
 * of the cfitsio library on loading a file into memory. Event selection is
 * kept optional (on compile time) as the data format may still evolve (for
 * the time being the TIME information is not valid).
 *
 * @todo Set data selection keywords.
 ***************************************************************************/
void ctselect::select(void)
{
    // Write Header for event selection
    if (logTerse()) {
        log << std::endl;
        log.header1("Event selection");
    }

    // Build selection string
    std::string selection;
    #if SELECT_TIME
    selection = "TIME >= "+str(m_tmin)+" && TIME <= "+str(m_tmax);
    if (logTerse())
        log << " Time range ................: " << m_tmin << "-" << m_tmax
            << std::endl;
    #endif
    #if SELECT_ENERGY
    if (selection.length() > 0)
        selection += " && ";
    selection += "ENERGY >= "+str(m_emin)+" && ENERGY <= "+str(m_emax);
    if (logTerse())
        log << " Energy range ..............: " << m_emin << "-" << m_emax
            << " TeV" << std::endl;
    #endif
    #if SELECT_ROI
    if (selection.length() > 0)
        selection += " && ";
    selection += "ANGSEP("+str(m_ra)+","+str(m_dec)+",RA,DEC) <= "+str(m_rad);
    if (logTerse()) {
        log << " Acceptance cone centre ....: RA=" << m_ra << ", DEC=" << m_dec
            << " deg" << std::endl;
        log << " Acceptance cone radius ....: " << m_rad << " deg" << std::endl;
    }
    #endif
    if (logTerse())
        log << " cfitsio selection .........: " << selection << std::endl;

    // Build input filename including selection expression
    std::string expression = m_infile;
    if (selection.length() > 0)
        expression += "[EVENTS]["+selection+"]";
    if (logTerse())
        log << " FITS filename .............: " << expression << std::endl;

    // Open FITS file
    m_file.open(expression);

    // Log selected FITS file
    if (logExplicit()) {
        log << std::endl;
        log.header1("FITS file content after selection");
        log << m_file << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append GTI to FITS file
 *
 * @todo We implement here a dummy method that overwrites an existing GTI
 *       entry but that only works if a single time interval exists. This
 *       is not clean, but since we have no method so far that allows
 *       deleting a HDU from a FITS file we have no other choice. Ideally,
 *       we would like to first delete an existing HDU and then append
 *       a new one. If this is not possible, we would like to delete the
 *       rows from an existing GTI and then fill it up anew. But in any
 *       case, the actual data selection method does not check the
 *       existing GTI, which we have to do to merge any subselection
 *       information in. In fact, we need a GTI merge method as well ...
 ***************************************************************************/
void ctselect::append_gti(void)
{
    // Setup selected time range
    GTime tstart;
    GTime tstop;
    tstart.met(m_tmin);
    tstop.met(m_tmax);

    // Initialise empty GTI table HDU
    GFitsTable* hdu = NULL;
    
    // Write selected time interval, either in an existing GTI or in a new
    // GTI that is appended to the FITS file
    try {
        hdu = m_file.table("GTI");
        GFitsTableDoubleCol* start = (GFitsTableDoubleCol*)hdu->column("START");
        GFitsTableDoubleCol* stop  = (GFitsTableDoubleCol*)hdu->column("STOP");
        (*start)(0) = tstart.met();
        (*stop)(0)  = tstop.met();
    }
    catch (GException::fits_hdu_not_found) {
        // Allocate Gti
        GGti gti;

        // Setup single GTI covering the selected time range
        gti.add(tstart, tstop);

        // Write GTI
        gti.write(&m_file);
    }

    // Log final FITS file
    if (logExplicit()) {
        log << std::endl;
        log.header1("Final FITS file content");
        log << m_file << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write data selection keywords
 *
 * @todo This is a very dumb data selection keyword writing routine that does
 *       not take into account any existing keywords. We definitely want a
 *       more secure logic that checks for existing keywords and possible
 *       conflicts. But for this prototype software, this code does the job.
 ***************************************************************************/
void ctselect::write_ds_keys(void)
{
    // Get event list header
    GFitsHDU* hdu = m_file.hdu("EVENTS");

    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Set cone selection string
        std::string dsval2 = "CIRCLE("+str(m_ra)+","+str(m_dec)+","+str(m_rad)+")";

        // Set energy selection string
        std::string dsval3 = str(m_emin)+":"+str(m_emax);

        // Add time selection keywords
        hdu->card("DSTYP1", "TIME",  "Data selection type");
        hdu->card("DSUNI1", "s",     "Data selection unit");
        hdu->card("DSVAL1", "TABLE", "Data selection value");
        hdu->card("DSREF1", ":GTI",  "Data selection reference");

        // Add acceptance cone selection
        hdu->card("DSTYP2", "POS(RA,DEC)", "Data selection type");
        hdu->card("DSUNI2", "deg",         "Data selection unit");
        hdu->card("DSVAL2", dsval2,        "Data selection value");
        
        // Add energy range selection
        hdu->card("DSTYP3", "ENERGY", "Data selection type");
        hdu->card("DSUNI3", "TeV",    "Data selection unit");
        hdu->card("DSVAL3", dsval3,   "Data selection value");

        // Set number of data selection keys
        hdu->card("NDSKEYS", 3,  "Number of data selections");

    } // endif: HDU was valid

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
void ctselect::init_members(void)
{
    // Initialise members
    m_infile.clear();
    m_outfile.clear();
    m_ra   = 0.0;
    m_dec  = 0.0;
    m_rad  = 0.0;
    m_tmin = 0.0;
    m_tmax = 0.0;
    m_emin = 0.0;
    m_emax = 0.0;

    // Set logger properties
    log.date(true);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctselect::free_members(void)
{
    // Free members

    // Return
    return;
}
