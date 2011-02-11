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
#include <cstdio>
#include "ctselect.hpp"
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
 * @brief Observations constructor
 *
 * This constructor creates an instance of the class that is initialised from
 * an observation container.
 ***************************************************************************/
ctselect::ctselect(GObservations obs) : GApplication(CTSELECT_NAME, CTSELECT_VERSION)
{
    // Initialise members
    init_members();

    // Set observations
    m_obs = obs;

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
 * @brief Copy constructor
 *
 * @param[in] app Application.
 ***************************************************************************/
ctselect::ctselect(const ctselect& app) : GApplication(app)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(app);

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

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] app Application.
 ***************************************************************************/
ctselect& ctselect::operator= (const ctselect& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->GApplication::operator=(app);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(app);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void ctselect::clear(void)
{
    // Free members
    free_members();
    this->GApplication::free_members();

    // Initialise members
    this->GApplication::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Execute application
 *
 * This method performs the event selection and saves the result
 ***************************************************************************/
void ctselect::execute(void)
{
    // Read ahead output filename so that it gets dumped correctly in the
    // parameters log
    m_outfile = (*this)["outfile"].filename();

    // Perform event selection
    run();

    // Save results
    save();

    // Return
    return;
}



/***********************************************************************//**
 * @brief Select event data
 *
 * This method reads in the application parameters and loops over all
 * observations that were found to perform an event selection. Event
 * selection is done by writing each observation to a temporary file and
 * re-opening the temporary file using the cfitsio event filter syntax.
 * The temporary file is deleted after this action so that no disk overflow
 * will occur. 
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

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        log.header1("Observations before selection");
        log << m_obs << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Event selection");
    }

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(&m_obs[i]);

        // Continue only if observation is a CTA observation
        if (obs != NULL) {

            // Write header for observation
            if (logTerse()) {
                if (obs->name().length() > 1)
                    log.header3("Observation "+obs->name());
                else
                    log.header3("Observation");
            }

            // Get temporary file name
            std::string filename = std::tmpnam(NULL);

            // Save observation in temporary file
            obs->save(filename, true);
            
            // Load observation from temporary file, including event selection
            select_events(obs, filename);

            // Remove temporary file
            std::remove(filename.c_str());

        } // endif: had a CTA observation

    } // endfor: looped over all observations

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        log.header1("Observations after selection");
        log << m_obs << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save simulated observation
 *
 * This method saves the results.
 ***************************************************************************/
void ctselect::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save observations");
    }

    // Get output filename
    m_outfile = (*this)["outfile"].value();

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(&m_obs[i]);

        // Save only if observation is a CTA observation
        if (obs != NULL) {

            // Save file
            obs->save(m_outfile, clobber());

        } // endif: observation was valid

    } // endfor: looped over observations

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
    // If there are no observations in container then add a single CTA
    // observation using the parameters from the parameter file
    if (m_obs.size() == 0) {

        // Get CTA event list file name
        m_infile = (*this)["infile"].filename();

        // Allocate CTA observation
        GCTAObservation obs;

        // Load CTA observation from file
        obs.load_unbinned(m_infile);

        // Append CTA observation to container
        m_obs.append(obs);

    } // endif: there was no observation in the container

    // Get parameters
    m_ra   = (*this)["ra"].real();
    m_dec  = (*this)["dec"].real();
    m_rad  = (*this)["rad"].real();
    m_tmin = (*this)["tmin"].real();
    m_tmax = (*this)["tmax"].real();
    m_emin = (*this)["emin"].real();
    m_emax = (*this)["emax"].real();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Select events
 *
 * Select events from a FITS file by making use of the selection possibility
 * of the cfitsio library on loading a file.
 *
 * @todo Implement a observation read() method to avoid saving the
 *       observation to a temporary file.
 ***************************************************************************/
void ctselect::select_events(GCTAObservation* obs, const std::string& filename)
{
    // Allocate selection string
    std::string selection;

    // Make time selection
    selection = "TIME >= "+str(m_tmin)+" && TIME <= "+str(m_tmax);
    if (logTerse())
        log << " Time range ................: " << m_tmin << "-" << m_tmax
            << std::endl;
    if (selection.length() > 0)
        selection += " && ";

    // Make energy selection
    selection += "ENERGY >= "+str(m_emin)+" && ENERGY <= "+str(m_emax);
    if (logTerse())
        log << " Energy range ..............: " << m_emin << "-" << m_emax
            << " TeV" << std::endl;
    if (selection.length() > 0)
        selection += " && ";

    // Make ROI selection
    selection += "ANGSEP("+str(m_ra)+","+str(m_dec)+",RA,DEC) <= "+str(m_rad);
    if (logTerse()) {
        log << " Acceptance cone centre ....: RA=" << m_ra << ", DEC=" << m_dec
            << " deg" << std::endl;
        log << " Acceptance cone radius ....: " << m_rad << " deg" << std::endl;
    }
    if (logTerse())
        log << " cfitsio selection .........: " << selection << std::endl;

    // Build input filename including selection expression
    std::string expression = filename;
    if (selection.length() > 0)
        expression += "[EVENTS]["+selection+"]";
    if (logTerse())
        log << " FITS filename .............: " << expression << std::endl;

    // Open FITS file
    GFits file(expression);

    // Log selected FITS file
    if (logExplicit()) {
        log << std::endl;
        log.header1("FITS file content after selection");
        log << file << std::endl;
    }

    // Get temporary file name
    std::string tmpname = std::tmpnam(NULL);
    
    // Save FITS file to temporary file
    file.saveto(tmpname, true);
            
    // Load observation from temporary file
    obs->load_unbinned(tmpname);

    // Get CTA event list pointer
    GCTAEventList* list = (GCTAEventList*)(obs->events());

    // Set ROI
    GCTARoi     roi;
    GCTAInstDir instdir;
    instdir.radec_deg(m_ra, m_dec);
    roi.centre(instdir);
    roi.radius(m_rad);
    list->roi(roi);

    // Set GTI
    GGti  gti;
    GTime tstart;
    GTime tstop;
    tstart.met(m_tmin);
    tstop.met(m_tmax);
    gti.append(tstart, tstop);
    list->gti(gti);

    // Set energy boundaries
    GEbounds ebounds;
    GEnergy  emin;
    GEnergy  emax;
    emin.TeV(m_emin);
    emax.TeV(m_emax);
    ebounds.append(emin, emax);
    list->ebounds(ebounds);

    // Remove temporary file
    std::remove(tmpname.c_str());

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
    m_obs.clear();
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
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctselect::copy_members(const ctselect& app)
{
    // Copy attributes
    m_infile  = app.m_infile;
    m_outfile = app.m_outfile;
    m_obs     = app.m_obs;
    m_ra      = app.m_ra;
    m_dec     = app.m_dec;
    m_rad     = app.m_rad;
    m_tmin    = app.m_tmin;
    m_tmax    = app.m_tmax;
    m_emin    = app.m_emin;
    m_emax    = app.m_emax;

    // Return
    return;
}    


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctselect::free_members(void)
{
    // Write separator into logger
    if (logTerse())
        log << std::endl;

    // Return
    return;
}
