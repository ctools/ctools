/***************************************************************************
 *                    ctselect - CTA data selection tool                   *
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
#define SELECT_TIME   0
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
 ***************************************************************************/
void ctselect::append_gti(void)
{
    // Allocate Gti
    GGti gti;

    // Setup single GTI covering the selected time range
    GTime tstart;
    GTime tstop;
    tstart.met(m_tmin);
    tstop.met(m_tmax);
    gti.add(tstart, tstop);

    // Write GTI
    gti.write(&m_file);

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
 * @brief Copy events
 *
 * This is some old code that is not used but that is kept for the moment
 * for debugging purposes.
 ***************************************************************************/
void ctselect::copy(void)
{
    // Open FITS file
    GFits file(m_infile);
    //std::cout << file << std::endl;

    // Optional perform Bit column testing
    #if BIT_COLUMN_TEST
    GFitsTableBitCol* ptr = (GFitsTableBitCol*)file.table("EVENTS")->column("TELMASK");
    for (int i = 0; i < ptr->length(); ++i) {
        std::cout << i << ": ";
        for (int k = 0; k < ptr->number(); ++k)
            std::cout << (int)(*(ptr))(i, k) << " ";
        std::cout << std::endl;
    }
    (*(ptr))(8507, 1) = true;
    (*(ptr))(8508, 0) = false;
    (*(ptr))(8509, 1) = false;
    (*(ptr))(8509, 2) = false;
    std::cout << "change now:" << std::endl;
    for (int i = 8507; i < 8510; ++i) {
        std::cout << i << ": ";
        for (int k = 0; k < ptr->number(); ++k)
            std::cout << (int)(*(ptr))(i, k) << " ";
        std::cout << std::endl;
    }
    #endif 

    // Optional perform String column testing
    #if STRING_COLUMN_TEST
    GFitsTableStringCol* ptr = (GFitsTableStringCol*)file.table("ANALYSIS")->column("UNIT");
    std::cout << "Anynul: " << ptr->anynul() << std::endl;
    for (int i = 0; i < ptr->length(); ++i) {
        std::cout << "\"" << ptr->string(i) << "\"" << std::endl;
    }
    #endif

    // Save file
    file.saveto(m_outfile, true);
    //std::cout << file << std::endl;

    // Optional perform String column testing
    #if STRING_COLUMN_TEST
    GFits file2(m_outfile);
    GFitsTableStringCol* ptr2 = (GFitsTableStringCol*)file2.table("ANALYSIS")->column("UNIT");
    std::cout << "Anynul: " << ptr2->anynul() << std::endl;
    for (int i = 0; i < ptr2->length(); ++i) {
        std::cout << "\"" << ptr2->string(i) << "\"" << std::endl;
    }
    #endif

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
