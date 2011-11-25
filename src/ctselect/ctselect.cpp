/***************************************************************************
 *                    ctselect - CTA data selection tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file ctselect.cpp
 * @brief CTA data selection tool implementation
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctselect.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_RUN                                               "ctselect::run()"

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

        // Initialise event filename
        m_infiles.push_back("");

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(&m_obs[i]);

        // Continue only if observation is a CTA observation
        if (obs != NULL) {

            // Write header for observation
            if (logTerse()) {
                if (obs->name().length() > 1) {
                    log.header3("Observation "+obs->name());
                }
                else {
                    log.header3("Observation");
                }
            }

            // Save event file name (for possible saving)
            m_infiles[i] = obs->eventfile();

            // Get temporary file name
            std::string filename = std::tmpnam(NULL);

            // Save observation in temporary file
            obs->save(filename, true);

            // Check temporary file
            std::string message = check_infile(filename);
            if (message.length() > 0) {
                throw GException::app_error(G_RUN, message);
            }

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
 * @brief Save observation(s) into event file(s)
 *
 * This method saves the observation(s) into event file(s). Each input event
 * file will be saved in a separate output event file. If a single input
 * event file is present, the "outfile" parameter specifies the filename of
 * the output event file. For multiple input event files, "outfile" is used
 * as suffix that is added to the input event file name to create the output
 * filename.
 ***************************************************************************/
void ctselect::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save observations");
    }

    // Get output filename (or prefix for multiple event files)
    m_outfile = (*this)["outfile"].value();

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(&m_obs[i]);

        // Save only if observation is a CTA observation
        if (obs != NULL) {

            // Get filename of events file
            std::string eventfile = m_infiles[i];

            // Set filename of output events file
            //TODO: we have to strip any path before appending the prefix
            //We may create a special GTools function for this.
            std::string outname = m_outfile;
            if (m_obs.size() > 1) {
                outname += "_" + eventfile;
                //TODO: strip .gz suffix
            }

            // Save observation into output events file
            obs->save(outname, clobber());

            // Copy over all other extensions that were present in the
            // input event file
            GFits infile(eventfile);
            GFits outfile(outname);
            for (int extno = 1; extno < infile.size(); ++extno) {
                GFitsHDU* hdu = infile.hdu(extno);
                if (hdu->extname() != "EVENTS" && hdu->extname() != "GTI") {
                    outfile.append(*hdu);
                }
            }
            infile.close();
            outfile.save(true);
            outfile.close();

        } // endif: observation was valid

    } // endfor: looped over observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. Times are assumed to be in the native CTA MJD format.
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
    m_expr = (*this)["expr"].string();

    // Derive time interval
    m_timemin.time(m_tmin, G_CTA_MJDREF, "days");
    m_timemax.time(m_tmax, G_CTA_MJDREF, "days");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Select events
 *
 * @param[in] obs CTA observation.
 * @param[in] filename File name.
 *
 * Select events from a FITS file by making use of the selection possibility
 * of the cfitsio library on loading a file. A selection string is created
 * from the specified criteria that is appended to the filename so that
 * cfitsio will automatically filter the event data. This selection string
 * is then applied when opening the FITS file. The opened FITS file is then
 * saved into a temporary file which is the loaded into the actual CTA
 * observation, overwriting the old CTA observation. The ROI, GTI and EBounds
 * of the CTA event list are then set accordingly to the specified selection.
 * Finally, the temporary file created during this process is removed.
 *
 * Good Time Intervals of the observation will be limited to the time
 * interval [m_tmin, m_tmax]. If m_tmin=m_tmax=0, no time selection is
 * performed.
 *
 * @todo Add an arbitrary selection string.
 * @todo Use INDEF instead of 0.0 for pointing as RA/DEC selection
 ***************************************************************************/
void ctselect::select_events(GCTAObservation* obs, const std::string& filename)
{
    // Allocate selection string
    std::string selection;
    char        cmin[80];
    char        cmax[80];
    char        cra[80];
    char        cdec[80];
    char        crad[80];

    // Set requested selections
    bool select_time = (m_tmin != 0.0 || m_tmax != 0.0);

    // Set RA/DEC selection
    double ra  = m_ra;
    double dec = m_dec;
    if (m_ra == 0.0 && m_dec == 0.0) {
        const GCTAPointing *pnt = obs->pointing();
        ra = pnt->dir().ra_deg();
        dec = pnt->dir().dec_deg();
    }

    // Set time selection interval. We make sure here that the time selection
    // interval cannot be wider than the GTIs covering the data. This is done
    // using GGti's reduce() method.
    if (select_time) {

        // Reduce GTIs to specified time interval. The complicated cast is
        // necessary here because the gti() method is declared const, so
        // we're not officially allowed to modify the GTIs.
        ((GGti*)(&obs->events()->gti()))->reduce(m_timemin, m_timemax);

    } // endif: time selection was required

    // Save GTI for later usage
    GGti gti = obs->events()->gti();

    // Make time selection
    if (select_time) {
    
        // Extract effective time interval in native MJD reference
        double tmin = gti.tstart().time();
        double tmax = gti.tstop().time();

        // Format time with sufficient accuracy and add to selection string
        sprintf(cmin, "%.8f", tmin);
        sprintf(cmax, "%.8f", tmax);
        selection = "TIME >= "+std::string(cmin)+" && TIME <= "+std::string(cmax);
        if (logTerse()) {
            log << " Time range ................: " << tmin << "-" << tmax
                << std::endl;
        }
        if (selection.length() > 0) {
            selection += " && ";
        }
    }

    // Make energy selection
    sprintf(cmin, "%.8f", m_emin);
    sprintf(cmax, "%.8f", m_emax);
    selection += "ENERGY >= "+std::string(cmin)+" && ENERGY <= "+std::string(cmax);
    if (logTerse()) {
        log << " Energy range ..............: " << m_emin << "-" << m_emax
            << " TeV" << std::endl;
    }
    if (selection.length() > 0) {
        selection += " && ";
    }

    // Make ROI selection
    sprintf(cra,  "%.6f", ra);
    sprintf(cdec, "%.6f", dec);
    sprintf(crad, "%.6f", m_rad);
    selection += "ANGSEP("+std::string(cra)+"," +
                 std::string(cdec)+",RA,DEC) <= " +
                 std::string(crad);
    if (logTerse()) {
        log << " Acceptance cone centre ....: RA=" << ra << ", DEC=" << dec
            << " deg" << std::endl;
        log << " Acceptance cone radius ....: " << m_rad << " deg" << std::endl;
    }
    if (logTerse())
        log << " cfitsio selection .........: " << selection << std::endl;

    // Add additional expression
    if (strip_whitespace(m_expr).length() > 0) {
        if (selection.length() > 0) {
            selection += " && ";
        }
        selection += "("+strip_whitespace(m_expr)+")";
    }

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
    instdir.radec_deg(ra, dec);
    roi.centre(instdir);
    roi.radius(m_rad);
    list->roi(roi);

    // Set GTI
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
    m_infiles.clear();
    m_timemin.clear();
    m_timemax.clear();
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
    m_infiles = app.m_infiles;
    m_timemin = app.m_timemin;
    m_timemax = app.m_timemax;
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


/***********************************************************************//**
 * @brief Check input filename
 *
 * @param[in] filename File name.
 *
 * This method checks if the input FITS file is correct.
 ***************************************************************************/
std::string ctselect::check_infile(const std::string& filename) const
{
    // Initialise message string
    std::string message = "";

    // Open FITS file
    GFits file(filename);

    // Check for EVENTS HDU
    GFitsTable* table = NULL;
    try {

        // Get pointer to FITS table
        table = file.table("EVENTS");

        // Initialise list of missing columns
        std::vector<std::string> missing;

        // Check for existence of TIME column
        if (!table->hascolumn("TIME")) {
            missing.push_back("TIME");
        }

        // Check for existence of ENERGY column
        if (!table->hascolumn("ENERGY")) {
            missing.push_back("ENERGY");
        }

        // Check for existence of RA column
        if (!table->hascolumn("RA")) {
            missing.push_back("RA");
        }

        // Check for existence of DEC column
        if (!table->hascolumn("DEC")) {
            missing.push_back("DEC");
        }

        // Set error message for missing columns
        if (missing.size() > 0) {
            message = "The following columns are missing in the"
                      " \"EVENTS\" extension of input file \""
                    + m_outfile + "\": ";
            for (int i = 0; i < missing.size(); ++i) {
                message += "\"" + missing[i] + "\"";
                if (i < missing.size()-1)
                    message += ", ";
            }
        }

    }
    catch (GException::fits_hdu_not_found& e) {
        message = "No \"EVENTS\" extension found in input file \""
                + m_outfile + "\".";
    }

    // Return
    return message;
}

