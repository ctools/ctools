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

/* __ Debug definitions __________________________________________________ */
#define BIT_COLUMN_TEST    0
#define STRING_COLUMN_TEST 1


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
 * @todo An interesting alternative is the fits_select_rows routine which copies
 * the selected rows from one FITS table to another. That's in fact the thing
 * we want to achieve here. The only way to use this routine would be to
 * open 2 FITS files, one for the input data and a fresh one for the output
 * data. The events table would then be copied over using the select routine,
 * any other tables would be copied over directly. The code would look like:
 *
 * GFits infile(m_infile);
 * GFits outfile(m_outfile);
 * for (int extno = 0; extno < infile.size(); ++extno) {
 *     GFitsHDU* inhdu  = infile.hdu(i);
 *     if (inhdu->extname() == "EVENTS") {
 *         GFitsHDU* outhdu = inhdu->select("expression");
 *         outfile.append_hdu(*outhdu);
 *     }
 *     else
 *         outfile.append_hdu(*inhdu);
 * }
 * outfile.save();
 * outfile.close();
 *
 ***************************************************************************/
void ctselect::select(void)
{
    // Open FITS file
    GFits file(m_infile);
    //std::cout << file << std::endl;

    // Optional perform Bit column testing
    #if BIT_COLUMN_TEST
    GFitsTableBitCol* ptr = (GFitsTableBitCol*)file.hdu("EVENTS")->column("TELMASK");
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
    GFitsTableStringCol* ptr = (GFitsTableStringCol*)file.hdu("ANALYSIS")->column("UNIT");
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
    GFitsTableStringCol* ptr2 = (GFitsTableStringCol*)file2.hdu("ANALYSIS")->column("UNIT");
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
