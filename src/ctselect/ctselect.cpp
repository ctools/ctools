/***************************************************************************
 *                      ctselect - Data selection tool                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2017 by Juergen Knoedlseder                         *
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
 * @brief Data selection tool definition
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>
#include "ctselect.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_RUN                                               "ctselect::run()"
#define G_SELECT_EVENTS          "ctselect::select_events(GCTAObservation*, "\
                                  "std::string&, std::string&, std::string&)"
#define G_SET_EBOUNDS    "ctselect::set_ebounds(GCTAObservation*, GEbounds&)"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_USE_MKSTEMP               //!< Use mkstemp for temporary filename


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
ctselect::ctselect(void) : ctobservation(CTSELECT_NAME, CTSELECT_VERSION)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observations constructor
 *
 * param[in] obs Observation container.
 *
 * Creates an instance of the class that is initialised using the information
 * provided in an observation container.
 ***************************************************************************/
ctselect::ctselect(const GObservations& obs) :
          ctobservation(CTSELECT_NAME, CTSELECT_VERSION, obs)
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
ctselect::ctselect(int argc, char *argv[]) : 
          ctobservation(CTSELECT_NAME, CTSELECT_VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Application.
 ***************************************************************************/
ctselect::ctselect(const ctselect& app) : ctobservation(app)
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
 * @return Application.
 ***************************************************************************/
ctselect& ctselect::operator=(const ctselect& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->ctobservation::operator=(app);

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
 * @brief Clear ctselect tool
 *
 * Clears ctselect tool.
 ***************************************************************************/
void ctselect::clear(void)
{
    // Free members
    free_members();
    this->ctobservation::free_members();
    this->ctool::free_members();

    // Clear base class (needed to conserve tool name and version)
    this->GApplication::clear();

    // Initialise members
    this->ctool::init_members();
    this->ctobservation::init_members();
    init_members();

    // Write header into logger
    log_header();

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
    if (logDebug()) {
        log.cout(true);
    }

    // Get parameters
    get_parameters();

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write header into logger
    log_header1(TERSE, "Event selection");

    // Initialise counters
    int n_observations = 0;

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Write header for the current observation
        log_header3(TERSE, get_obs_header(m_obs[i]));

        // Initialise event input and output filenames and the event
        // and GTI extension names
        m_infiles.push_back("");
        m_evtname.push_back(gammalib::extname_cta_events);
        m_gtiname.push_back(gammalib::extname_gti);

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Skip observation if it's not CTA
        if (obs == NULL) {
            std::string msg = " Skipping "+m_obs[i]->instrument()+
                              " observation";
            log_string(NORMAL, msg);
            continue;
        }

        // Skip observation if we have a binned observation
        if (obs->eventtype() == "CountsCube") {
            std::string msg = " Skipping binned "+m_obs[i]->instrument()+
                              " observation";
            log_string(NORMAL, msg);
            continue;
        }

        // Increment counter
        n_observations++;

        // Save event file name (for possible saving)
        m_infiles[i] = obs->eventfile();

        // Extract event and GTI extension names from input FITS file
        GFilename fname(m_infiles[i]);
        if (fname.has_extname()) {
            m_evtname[i] = fname.extname();
        }
        m_gtiname[i] = get_gtiname(fname.url(), m_evtname[i]);

        // Write input file information into logger
        log_value(NORMAL, "Input filename", m_infiles[i]);
        log_value(NORMAL, "Event extension name", m_evtname[i]);
        log_value(NORMAL, "GTI extension name", m_gtiname[i]);

        // Fall through in case that the event file is empty
        if (obs->events()->size() == 0) {
            log_string(NORMAL, " Warning: No events in event file \""+
                       m_infiles[i]+"\". Event selection skipped.");
            continue;
        }

        // If we have an input file then check it
        if (!m_infiles[i].empty()) {
            std::string message = check_infile(m_infiles[i], m_evtname[i]);
            if (!message.empty()) {
                throw GException::invalid_value(G_RUN, message);
            }
        }

        // Get temporary file name
        #if G_USE_MKSTEMP
        char tpl[]  = "ctselectXXXXXX";
        int  fileid = mkstemp(tpl);
        std::string filename(tpl);
        #else
        std::string filename = std::tmpnam(NULL);
        #endif

        // Save observation in temporary file. We add here the events and
        // GTI extension name so that the GCTAObservation::save method can
        // use this information for writing the proper extension names into
        // the temporary file
        obs->save(filename+"["+m_evtname[i]+";"+ m_gtiname[i]+"]", true);

        // Log saved FITS file.
        if (logExplicit()) {
            GFits tmpfile(filename);
            log.header3("FITS file content of temporary file");
            log << tmpfile << std::endl;
            tmpfile.close();
        }

        // If we have a temporary file then check it
        if (!filename.empty()) {
            std::string message = check_infile(filename, m_evtname[i]);
            if (!message.empty()) {
                throw GException::invalid_value(G_RUN, message);
            }
        }

        // Load observation from temporary file, including event selection
        select_events(obs, filename, m_evtname[i], m_gtiname[i]);

        // Close temporary file
        #if G_USE_MKSTEMP
        close(fileid);
        #endif

        // Remove temporary file
        std::remove(filename.c_str());

    } // endfor: looped over all observations

    // If more than a single observation has been handled then make sure that
    // an XML file will be used for storage
    if (n_observations > 1) {
        m_use_xml = true;
    }

    // Write observation(s) into logger
    log_observations(NORMAL, m_obs, "Output observation");

    // Optionally publish event list(s)
    if ((*this)["publish"].boolean()) {
        publish();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save the selected event list(s)
 *
 * This method saves the selected event list(s) into FITS file(s). There are
 * two modes, depending on the m_use_xml flag.
 *
 * If m_use_xml is true, all selected event list(s) will be saved into FITS
 * files, where the output filenames are constructued from the input
 * filenames by prepending the m_prefix string to name. Any path information
 * will be stripped form the input name, hence event files will be written
 * into the local working directory (unless some path information is present
 * in the prefix). In addition, an XML file will be created that gathers
 * the filename information for the selected event list(s). If an XML file
 * was present on input, all metadata information will be copied from this
 * input file.
 *
 * If m_use_xml is false, the selected event list will be saved into a FITS
 * file.
 ***************************************************************************/
void ctselect::save(void)
{
    // Write header into logger
    log_header1(TERSE, gammalib::number("Save event list", m_obs.size()));

    // Case A: Save event file(s) and XML metadata information
    if (m_use_xml) {
        save_xml();
    }

    // Case B: Save event file as FITS file
    else {
        save_fits();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Publish event lists
 *
 * @param[in] name Event list name.
 ***************************************************************************/
void ctselect::publish(const std::string& name)
{
    // Write header into logger
    log_header1(TERSE, gammalib::number("Publish event list", m_obs.size()));

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Handle only CTA observations
        if (obs != NULL) {

            // Continue only if there is an event list
            if (obs->events()->size() != 0) {

                // Set default name if user name is empty
                std::string user_name(name);
                if (user_name.empty()) {
                    user_name = CTSELECT_NAME;
                }

                // If there are several event lists then add an index
                if (m_use_xml) {
                    user_name += gammalib::str(i);
                }

                // Write event list name into logger
                log_value(NORMAL, "Event list name", user_name);

                // Write events into in-memory FITS file
                GFits fits;
                obs->write(fits);

                // Publish
                fits.publish(gammalib::extname_cta_events, user_name);

            } // endif: there were events

        } // endif: observation was a CTA observation

    } // endfor: looped over observations

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
    // Initialise parameters
    m_outobs.clear();
    m_prefix.clear();
    m_usepnt = false;
    m_roi.clear();
    m_tmin   = 0.0;
    m_tmax   = 0.0;
    m_emin   = 0.0;
    m_emax   = 0.0;
    m_expr.clear();
    m_usethres.clear();
    m_chatter = static_cast<GChatter>(2);

    // Initialise protected members
    m_infiles.clear();
    m_evtname.clear();
    m_gtiname.clear();
    m_timemin.clear();
    m_timemax.clear();
    m_select_energy = true;
    m_select_roi    = true;
    m_select_time   = true;

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
    // Copy parameters
    m_outobs   = app.m_outobs;
    m_prefix   = app.m_prefix;
    m_usepnt   = app.m_usepnt;
    m_roi      = app.m_roi;
    m_tmin     = app.m_tmin;
    m_tmax     = app.m_tmax;
    m_emin     = app.m_emin;
    m_emax     = app.m_emax;
    m_expr     = app.m_expr;
    m_usethres = app.m_usethres;
    m_chatter  = app.m_chatter;

    // Copy protected members
    m_infiles       = app.m_infiles;
    m_evtname       = app.m_evtname;
    m_gtiname       = app.m_gtiname;
    m_timemin       = app.m_timemin;
    m_timemax       = app.m_timemax;
    m_select_energy = app.m_select_energy;
    m_select_roi    = app.m_select_roi;
    m_select_time   = app.m_select_time;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctselect::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all user parameters from parameter file or (if required) by querying
 * the user. Times are assumed to be in the native CTA MJD format.
 *
 * This method also loads observations if no observations are yet allocated.
 * Observations are either loaded from a single CTA even list, or from a
 * XML file using the metadata information that is stored in that file.
 ***************************************************************************/
void ctselect::get_parameters(void)
{
    // Initialise selection flags
    m_select_energy = true;
    m_select_roi    = true;
    m_select_time   = true;

    // Setup observations from "inobs" parameter. Do not request response
    // information and do not accept counts cubes.
    setup_observations(m_obs, false, true, false);

    // Get some parameters
    m_usepnt = (*this)["usepnt"].boolean();

    // Get the RoI and enable RoI selection if the RoI is valid
    m_roi        = get_roi();
    m_select_roi = m_roi.is_valid();

    // Check for sanity of time selection parameters
    if ((*this)["tmin"].is_valid() && (*this)["tmax"].is_valid()) {

        // Get User parameters
        m_tmin = (*this)["tmin"].real();
        m_tmax = (*this)["tmax"].real();

        // Additional check for time values
        if (m_tmin >= m_tmax) {
            m_select_time = false;
        }
        else {
            m_select_time = true;
        }
    }
    else {
        m_select_time = false;
    }

    // Check for sanity of energy selection parameters
    if ((*this)["emin"].is_valid() && (*this)["emax"].is_valid()) {
        m_emin          = (*this)["emin"].real();
        m_emax          = (*this)["emax"].real();
        m_select_energy = true;
    }
    else {
        m_select_energy = false;
    }

    // Get other User parameters
    m_expr     = (*this)["expr"].string();
    m_usethres = (*this)["usethres"].string();
    m_chatter  = static_cast<GChatter>((*this)["chatter"].integer());

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outobs = (*this)["outobs"].filename();
        m_prefix = (*this)["prefix"].string();
    }

    // Set time interval with input times given in CTA reference
    // time (in seconds)
    m_timemin.set(m_tmin, m_cta_ref);
    m_timemax.set(m_tmax, m_cta_ref);

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Select events
 *
 * @param[in,out] obs CTA observation.
 * @param[in] filename File name.
 * @param[in] evtname Event extension name.
 * @param[in] gtiname GTI extension name.
 *
 * @exception GException::invalid_value
 *            No events extension found in FITS file.
 *
 * Select events from a FITS file by making use of the selection possibility
 * of the cfitsio library on loading a file. A selection string is created
 * from the specified criteria that is appended to the filename so that
 * cfitsio will automatically filter the event data. This selection string
 * is then applied when opening the FITS file. The event list in the current
 * observation is replaced by selected event list read from the FITS file.
 *
 * Good Time Intervals of the observation will be limited to the time
 * interval [m_tmin, m_tmax]. If m_tmin=m_tmax=0, no time selection is
 * performed.
 ***************************************************************************/
void ctselect::select_events(GCTAObservation*   obs,
                             const std::string& filename,
                             const std::string& evtname,
                             const std::string& gtiname)
{
    // Write header into logger
    log_header3(NORMAL, "Events selection");

    // Initialise selection and addition strings
    std::string selection;
    std::string add;

    // Initialise selection flags
    bool remove_all = false;

    // Get CTA event list pointer
    GCTAEventList* list =
        static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));

    // Get existing Roi and energy bounds for possible later use
    // (will be empty if unavailable)
    GCTARoi  old_roi     = list->roi();
    GEbounds old_ebounds = list->ebounds();

    // Determine new energy boundaries for selection
    // taking into account previous existing ones
    double   emin    = 0.0;
    double   emax    = 0.0;
    GEbounds ebounds = set_ebounds(obs, list->ebounds());
    if (ebounds.size() >= 1) {
        emin = ebounds.emin(0).TeV();
        emax = ebounds.emax(0).TeV();
    }

    // Analyse expression to see if a selection is required
    bool select_expr = (gammalib::strip_whitespace(m_expr).length() > 0);

    // Set RoI selection. If the "usepnt" parameter is set to true then
    // use the pointing direction as the RoI centre.
    double ra  = m_roi.centre().dir().ra_deg();
    double dec = m_roi.centre().dir().dec_deg();
    double rad = m_roi.radius();
    if (m_usepnt) {
        const GCTAPointing &pnt = obs->pointing();
        ra  = pnt.dir().ra_deg();
        dec = pnt.dir().dec_deg();
    }
    
    // Set time selection interval. We make sure here that the time selection
    // interval cannot be wider than the GTIs covering the data. This is done
    // using GGti's reduce() method.
    if (m_select_time) {

        // Reduce GTIs to specified time interval. The complicated cast is
        // necessary here because the gti() method is declared const, so
        // we're not officially allowed to modify the GTIs.
        ((GGti*)(&list->gti()))->reduce(m_timemin, m_timemax);

    } // endif: time selection was required

    // Save GTI for later usage
    GGti gti = list->gti();

    // Make time selection
    if (m_select_time) {
    
        // Extract effective time interval in the reference time of the
        // event list. We get this reference time from gti.reference().
        double tmin = gti.tstart().convert(gti.reference());
        double tmax = gti.tstop().convert(gti.reference());

        // Format time with sufficient accuracy and add to selection string
        char cmin[80];
        char cmax[80];
        sprintf(cmin, "%.8f", tmin);
        sprintf(cmax, "%.8f", tmax);
        selection = "TIME >= "+std::string(cmin)+" && TIME <= "+std::string(cmax);
        add       = " && ";
        log_value(NORMAL, "Time range",
                  gammalib::str(tmin)+" - "+gammalib::str(tmax)+" s");

    } // endif: made time selection

    // Make energy selection
    if (m_select_energy) {

        // Log the requested energy selection
        if (logNormal()) {
            log << gammalib::parformat("Selected energy range");
            if (emin > 0.0 && emax > 0.0) {
                if (emin >= emax) {
                    log << "None. There is no overlap between existing ";
                    log << "and requested energy range." << std::endl;
                }
                else {
                    log << emin << " - " << emax << " TeV" << std::endl;
                }
            }
            else if (emin > 0.0) {
                log << "> " << emin << " TeV" << std::endl;
            }
            else if (emax > 0.0) {
                log << "< " << emax << " TeV" << std::endl;
            }
            else {
                log << "None" << std::endl;
            }
        }

        // Apply the requested energy selection
        if (emin > 0.0 && emax > 0.0) {
            if (emin >= emax) {
                remove_all = true;
            }
            else {
                char cmin[80];
                char cmax[80];
                sprintf(cmin, "%.8f", emin);
                sprintf(cmax, "%.8f", emax);
                selection += add + "ENERGY >= " + std::string(cmin)+
                                   " && ENERGY <= " + std::string(cmax);
                add        = " && ";
            }
        }
        else if (emin > 0.0) {
            char cmin[80];
            sprintf(cmin, "%.8f", emin);
            selection += add + "ENERGY >= " + std::string(cmin);
            add        = " && ";
        }
        else if (emax > 0.0) {
            char cmax[80];
            sprintf(cmax, "%.8f", emax);
            selection += add + "ENERGY <= " + std::string(cmax);
            add        = " && ";
        }
        else {
            remove_all = true;
        }

    } // endif: made energy selection

    // Make RoI selection
    if (m_select_roi) {

        // Store the original radius
        double original_rad = rad;

        // Log the requested RoI
        log_value(NORMAL, "Requested RoI",
                  "Centre(RA,DEC)=("+gammalib::str(ra)+", "+
                  gammalib::str(dec)+") deg, Radius="+gammalib::str(rad)+
                  " deg");

        // If we have already an RoI then make sure that the selected
        // ROI overlaps with the existing RoI
        double roi_radius = list->roi().radius();
        if (roi_radius > 0.0) {
            GSkyDir roi_centre = list->roi().centre().dir();
            log_value(NORMAL, "RoI of data",
                      "Centre(RA,DEC)=("+gammalib::str(roi_centre.ra_deg())+
                      ", "+gammalib::str(roi_centre.dec_deg())+") deg, Radius="+
                      gammalib::str(roi_radius)+" deg");
            GSkyDir centre;
            centre.radec_deg(ra, dec);
            double distance = centre.dist_deg(roi_centre);
            if (distance + rad > roi_radius) {
                rad = roi_radius - distance;
            }
        }

        // If the RoI radius is negative then there is no overlap between
        // the existing and the requested RoI and the radius will be restored
        // to the original value. We signal that all events should be
        // removed.
        if (rad <= 0.0) {
            rad        = original_rad;
            remove_all = true;
            log_value(NORMAL, "Selected RoI", "None. There is no overlap "
                      "between existing and requested RoI.");
        }
        else {
            log_value(NORMAL, "Selected RoI", "Centre(RA,DEC)=("+
                      gammalib::str(ra)+", "+gammalib::str(dec)+") deg, "
                      "Radius="+gammalib::str(rad)+" deg");
        }
        
        // Format RoI selection
        char cra[80];
        char cdec[80];
        char crad[80];
        sprintf(cra,  "%.6f", ra);
        sprintf(cdec, "%.6f", dec);
        sprintf(crad, "%.6f", rad);
        selection += add + "ANGSEP("+std::string(cra)+"," +
                           std::string(cdec)+",RA,DEC) <= " +
                           std::string(crad);
        add        = " && ";

    } // endif: made ROI selection

    // Make an expression selection
    if (select_expr) {

        // Append the expression
        selection += add + "("+gammalib::strip_whitespace(m_expr)+")";

    } // endif: made expression selection

    // Dump cfitsio selection string
    log_value(NORMAL, "cfitsio selection", selection);

    // Build input filename including selection expression
    std::string expression = filename + "[" + evtname + "]";
    if (selection.length() > 0) {
        expression += "["+selection+"]";
    }

    // Dump FITS filename including selection expression
    log_value(NORMAL, "FITS filename", expression);

    // Open FITS file
    GFits file(expression);

    // Log selected FITS file
    log_header3(EXPLICIT, "FITS file content after selection");
    log_string(EXPLICIT, file.print(m_chatter));

    // Check if we have an events HDU
    if (!file.contains(evtname)) {
        std::string msg = "No events extension \""+evtname+"\" found in "
                          "FITS file. The expression \""+expression+"\" "
                          "was used to open the FITS file.";
        throw GException::invalid_value(G_SELECT_EVENTS, msg);
    }

    // Determine number of events in the events HDU. If removal of all events
    // has been requested then set the number of events to zero.
    int nevents = (remove_all) ? 0 : file.table(evtname)->nrows();

    // If the selected event list is empty then append an empty event list
    // to the observation
    if (nevents < 1) {

        // Create empty event list
        GCTAEventList eventlist;

        // Append list to observation
        obs->events(eventlist);

    }

    // ... otherwise load the data from the temporary file
    else {
        obs->read(file);
    }

    // Get CTA event list pointer
    list = static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));

    // Make sure that events are fetched since the temporary file will be
    // closed later
    list->fetch();

    // If RoI selection has been applied then set the event list RoI
    if (m_select_roi) {
        GCTAInstDir instdir;
        instdir.dir().radec_deg(ra, dec);
        list->roi(GCTARoi(instdir, rad));
    }

    // ... otherwise restore old RoI information if it existed
    else if (old_roi.is_valid()) {
        list->roi(old_roi);
    }

    // Set event list GTI (in any case as any event list has a GTI)
    list->gti(gti);

    // If an energy selection has been applied then set the energy boundaries
    if (m_select_energy) {
        GEbounds ebounds;
        ebounds.append(GEnergy(emin, "TeV"), GEnergy(emax, "TeV"));
        list->ebounds(ebounds);
    }

    // ... otherwise restore old Ebounds if they existed
    else if (old_ebounds.size() > 0) {
        list->ebounds(old_ebounds);
    }

    // Recompute ontime and livetime.
    obs->ontime(list->gti().ontime());
    obs->livetime(list->gti().ontime() * obs->deadc());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return energy boundaries for a given observation
 *
 * @param[in] obs Pointer to CTA observation.
 * @param[in] ebounds Current energy boundaries.
 * @return Energy boundaries.
 *
 * Returns the energy boundaries for a given observation. Depending on the
 * value of the usethres parameter, the following values will be
 * returned:
 *
 *     NONE:    [max(emin,emin_exist),min(emax,emax_exist)]
 *     DEFAULT: [max(emin,emin_exist,emin_save),min(emax,emax_exist,emax_save)]
 *     USER:    [max(emin,emin_exist,emin_user),min(emax,emax_exist,emax_user)]
 *
 * where
 *
 *     emin is the value of the emin parameter
 *     emax is the value of the emax parameter
 *     emin_exist is the value of any existing minimum boundary
 *     emax_exist is the value of any existing maximum boundary
 *     emin_save is the lower save threshold
 *     emax_save is the upper save threshold
 *     emin_user is the lower user threshold
 *     emax_user is the upper user threshold
 *
 * Any threshold value of 0 will be ignored.
 ***************************************************************************/
GEbounds ctselect::set_ebounds(GCTAObservation* obs, const GEbounds& ebounds) const
{
    // Set emin and emax values from user parameters
    double emin = m_emin;
    double emax = m_emax;

    // Check if we have already energy boundaries
    if (ebounds.size() > 0) {

        // Raise minimum energy to lower boundary
        if ((emin == 0.0) ||
            ((ebounds.emin(0).TeV() > 0.0) && (emin < ebounds.emin(0).TeV()))) {
            emin = ebounds.emin(0).TeV();
        }

        // Lower maximum energy to upper boundary
        if ((emax == 0.0) ||
            ((ebounds.emax(0).TeV() > 0.0) && (emax > ebounds.emax(0).TeV()))) {
            emax = ebounds.emax(0).TeV();
        }

    } // endif: there were already energy boundaries

    // Check if default threshold (the one from the IRF, also known as save
    // threshold) should be applied
    if (m_usethres == "DEFAULT") {

        // Get CTA IRF respsonse pointer and throw exception if we don't
        //  find IRF
        const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>
                                                             (obs->response());
        if (rsp == NULL) {
            std::string msg = "No IRF response attached to given observation.";
            throw GException::invalid_value(G_SET_EBOUNDS, msg);
        }

        // Retrieve energy range from response information
        double lo_save_thres = rsp->lo_save_thres();
        double hi_save_thres = rsp->hi_save_thres();

        // Check threshold information for validity
        if ((lo_save_thres > 0.0) &&
            (hi_save_thres > 0.0) &&
            (lo_save_thres >= hi_save_thres)) {
            std::string msg = "IRF \""+obs->name()+"\" contains an invalid "
                              "energy range ["+gammalib::str(lo_save_thres)+","+
                              gammalib::str(hi_save_thres)+"] TeV.";
            throw GException::invalid_value(G_SET_EBOUNDS, msg);
        }

        // Raise minimum energy to lower threshold (if lower threshold exists)
        if ((emin == 0.0) ||
            ((lo_save_thres > 0.0) && (emin < lo_save_thres))) {
            emin = lo_save_thres;
        }

        // Lower maximum energy to upper threshold (if upper threshold exists)
        if ((emax == 0.0) ||
            ((hi_save_thres > 0.0) && (emax > hi_save_thres))) {
            emax = hi_save_thres;
        }

    } //endif: usethres was "DEFAULT"

    // ... otherwise check if a user specified threshold should be applied
    else if (m_usethres == "USER") {

        // Retrieve thresholds from observation
        double lo_user_thres = obs->lo_user_thres();
        double hi_user_thres = obs->hi_user_thres();

        // Check threshold information for validity
        if ((lo_user_thres > 0.0) &&
            (hi_user_thres > 0.0) &&
            (lo_user_thres >= hi_user_thres)) {
            std::string msg = "User energy range ["+gammalib::str(lo_user_thres)+
                              ","+gammalib::str(hi_user_thres)+"] TeV is invalid.";
            throw GException::invalid_value(G_SET_EBOUNDS, msg);
        }

        // Raise minimum energy to lower threshold (if lower threshold exists)
        if ((emin == 0.0) ||
            ((lo_user_thres > 0.0) && (emin < lo_user_thres))) {
            emin = lo_user_thres;
        }

        // Lower maximum energy to upper threshold (if upper threshold exists)
        if ((emax == 0.0) ||
            ((hi_user_thres > 0.0) && (emax > hi_user_thres))) {
            emax = hi_user_thres;
        }

    } // endif: usethres was "USER

    // Set selection energy boundaries
    GEbounds result;
    if (emax > emin) {
        result.append(GEnergy(emin, "TeV"), GEnergy(emax, "TeV"));
    }

    // Return result
    return result;

}


/***********************************************************************//**
 * @brief Check input filename
 *
 * @param[in] filename File name.
 * @param[in] evtname Event extension name.
 *
 * This method checks if the input FITS file is correct.
 ***************************************************************************/
std::string ctselect::check_infile(const std::string& filename,
                                   const std::string& evtname) const
{
    // Initialise message string
    std::string message = "";

    // Open FITS file
    GFits fits(filename);

    // Check for existence of events extensions
    if (!fits.contains(evtname)) {
        message = "No \""+evtname+"\" extension found in input file \""+
                  filename + "\".";
    }

    // ... otherwise check column names
    else {

        // Get pointer to FITS table
        GFitsTable* table = fits.table(evtname);

        // Initialise list of missing columns
        std::vector<std::string> missing;

        // Check for existence of TIME column
        if (!table->contains("TIME")) {
            missing.push_back("TIME");
        }

        // Check for existence of ENERGY column
        if (!table->contains("ENERGY")) {
            missing.push_back("ENERGY");
        }

        // Check for existence of RA column
        if (!table->contains("RA")) {
            missing.push_back("RA");
        }

        // Check for existence of DEC column
        if (!table->contains("DEC")) {
            missing.push_back("DEC");
        }

        // Set error message for missing columns
        if (!missing.empty()) {
            message = "The following columns are missing in the "
                      "\""+evtname+"\" extension of input file \""+
                      filename + "\": ";
            for (int i = 0; i < missing.size(); ++i) {
                message += "\"" + missing[i] + "\"";
                if (i < missing.size()-1) {
                    message += ", ";
                }
            }
        }

    } // endelse: checked column names

    // Return
    return message;
}


/***********************************************************************//**
 * @brief Set output file name.
 *
 * @param[in] filename Input file name.
 *
 * Converts an input file name into an output filename by prepending the
 * prefix stored in the member m_prefix to the input file name. Any path as
 * well as extension will be stripped from the input file name. Also a
 * trailing ".gz" will be stripped as one cannot write into gzipped files.
 ***************************************************************************/
std::string ctselect::set_outfile_name(const std::string& filename) const
{
    // Create filename
    GFilename fname(filename);

    // Split input filename without any extensions into path elements
    std::vector<std::string> elements = gammalib::split(fname.url(), "/");

    // The last path element is the filename
    std::string outname = m_prefix + elements[elements.size()-1];

    // Strip any ".gz"
    outname = gammalib::strip_chars(outname, ".gz");
    
    // Return output filename
    return outname;
}


/***********************************************************************//**
 * @brief Get Good Time Intervals extension name
 *
 * @param[in] filename Input file name.
 * @param[in] evtname Events extension name.
 *
 * Extracts the Good Time Intervals extension name from the event file. We
 * do this by loading the events and accessing the Good Time Intervals
 * extension name using the GCTAEventList::gtiname() method. If the file name
 * is empty, the method returns `GTI`.
 ***************************************************************************/
std::string ctselect::get_gtiname(const std::string& filename,
                                  const std::string& evtname) const
{
    // Initialise GTI name
    std::string gtiname = gammalib::extname_gti;

    // Continue only if the filename is not empty
    if (!filename.empty()) {

        // Load events
        GCTAEventList events(filename+"["+evtname+"]");

        // Get GTI name
        gtiname = events.gtiname();

    }

    // Return GTI name
    return (gtiname);
}


/***********************************************************************//**
 * @brief Save event list in FITS format.
 *
 * Save the event list as a FITS file. The file name of the FITS file is
 * specified by the "outobs" parameter.
 ***************************************************************************/
void ctselect::save_fits(void)
{
    // Save only if there are observations
    if (m_obs.size() > 0) {

        // Get output filename
        m_outobs = (*this)["outobs"].filename();

        // Get CTA observation from observation container
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

        // Save only if it's a CTA observation
        if (obs != NULL) {
    
            // Save only if file name is non-empty
            if (m_infiles[0].length() > 0) {

                // Create file name object
                GFilename fname(m_outobs);

                // Extract filename and event extension name
                std::string outfile = fname.url();

                // Append event extension name. We handle here the possibility
                // to write the events into a different extension.
                if (fname.has_extname()) {
                    outfile += "["+fname.extname()+"]";
                }
                else {
                    outfile += "["+m_evtname[0]+"]";
                }

                // Log filename
                log_value(NORMAL, "Event list file", outfile);

                // Save event list
                save_event_list(obs, m_infiles[0], m_evtname[0], m_gtiname[0],
                                outfile);

            } // endif: filename was non empty

        } // endif: observation was CTA observation

    } // endif: there were observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save event list(s) in XML format.
 *
 * Save the event list(s) into FITS files and write the file path information
 * into a XML file. The filename of the XML file is specified by the outfile
 * parameter, the filename(s) of the event lists are built by prepending a
 * prefix to the input event list filenames. Any path present in the input
 * filename will be stripped, i.e. the event list(s) will be written in the
 * local working directory (unless a path is specified in the prefix).
 ***************************************************************************/
void ctselect::save_xml(void)
{
    // Get output filename and prefix
    m_outobs = (*this)["outobs"].filename();
    m_prefix = (*this)["prefix"].string();

    // Issue warning if output filename has no .xml suffix
    log_string(TERSE, warn_xml_suffix(m_outobs));

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Skip observations that are no CTA observations
        if (obs == NULL) {
            continue;
        }

        // Skip observations that have empty names
        if (m_infiles[i].length() == 0) {
            continue;
        }

        // Set event output file name
        std::string outfile = set_outfile_name(m_infiles[i]);

        // Append event extension name
        outfile += "["+m_evtname[i]+"]";

        // Log filename
        log_value(NORMAL, "Event list file", outfile);

        // Store output file name in observation
        obs->eventfile(outfile);

        // Save event list
        save_event_list(obs, m_infiles[i], m_evtname[i], m_gtiname[i],
                        outfile);

    } // endfor: looped over observations

    // Write observation definition XML file name into logger
    log_value(NORMAL, "Obs. definition file", m_outobs);

    // Save observations in XML file
    m_obs.save(m_outobs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save event list into FITS file
 *
 * @param[in] obs Pointer to CTA observation.
 * @param[in] infile Input file name.
 * @param[in] evtname Event extension name.
 * @param[in] gtiname GTI extension name.
 * @param[in] outfile Output file name.
 *
 * Saves an event list and the corresponding Good Time Intervals into a FITS
 * file and copy all others extensions from the input file to the output
 * file.
 *
 * If an extension name is specified in the @p outfile argument, the events
 * and eventually also the Good Time Intervals will be extracted from the
 * argument and used for writing the events. The format is
 *
 *      <filename>[<event extension name;GTI extension name>]
 *
 * where <filename> needs to be replaced by the name of the FITS file,
 * and <event extension name;GTI extension name> by the name of the events
 * and Good Time Intervals extensions. For example 
 *
 *      myfits.fits[EVENTS1;GTI1]
 *
 * will write the selected events into the "EVENTS1" extension and the
 * Good Time Intervals into the "GTI1" extension of the "myfits.fits" FITS
 * file. If the Good Time Intervals extension name is skipped, e.g.
 *
 *      myfits.fits[EVENTS1]
 *
 * the original extension name for the Good Time Intervals will be kept.
 * Analogously, only the Good Time Intervals extension name can be changed
 * by specifying
 *
 *      myfits.fits[;GTI1]
 *
 * In none of the cases will the original events and Good Time Intervals be
 * copied over to the output file.
 ***************************************************************************/
void ctselect::save_event_list(const GCTAObservation* obs,
                               const std::string&     infile,
                               const std::string&     evtname,
                               const std::string&     gtiname,
                               const std::string&     outfile) const
{
    // Save only if we have an event list
    if (obs->eventtype() == "EventList") {

        // Set output FITS file event extension names
        GFilename   outname(outfile);
        std::string outevt = evtname;
        std::string outgti = gtiname;
        if (outname.has_extname()) {
            std::vector<std::string> extnames =
                       gammalib::split(outname.extname(), ";");
            if (extnames.size() > 0) {
                std::string extname = gammalib::strip_whitespace(extnames[0]);
                if (!extname.empty()) {
                    outevt = extname;
                }
            }
            if (extnames.size() > 1) {
                std::string extname = gammalib::strip_whitespace(extnames[1]);
                if (!extname.empty()) {
                    outgti = extname;
                }
            }
        }

        // Create output FITS file
        GFits outfits;

        // Write observation into FITS file
        obs->write(outfits, outevt, outgti);

        // Copy all extensions other than evtname and gtiname extensions
        // from the input to the output event list. The evtname and
        // gtiname extensions are written by the save method, all others
        // that may eventually be present have to be copied over
        // explicitly.
        GFits infits(infile);
        for (int extno = 1; extno < infits.size(); ++extno) {
            GFitsHDU* hdu = infits.at(extno);
            if (hdu->extname() != evtname &&
                hdu->extname() != gtiname &&
                hdu->extname() != outevt  &&
                hdu->extname() != outgti) {
                outfits.append(*hdu);
            }
        }

        // Close input file
        infits.close();

        // Save file to disk and close it (we need both operations)
        outfits.saveto(outname.url(), clobber());
        outfits.close();

    } // endif: observation was unbinned

    // Return
    return;
}
