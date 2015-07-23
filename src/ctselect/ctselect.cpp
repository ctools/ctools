/***************************************************************************
 *                      ctselect - Data selection tool                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2015 by Juergen Knoedlseder                         *
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
                                                              "std::string&)"
#define G_SET_EBOUNDS    "ctselect::set_ebounds(GCTAObservation*, GEbounds&)"
#define G_GET_PARAMETERS                         "ctselect::get_parameters()"

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
ctselect::ctselect(void) : ctool(CTSELECT_NAME, CTSELECT_VERSION)
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
 * This constructor creates an instance of the class that is initialised from
 * an observation container.
 ***************************************************************************/
ctselect::ctselect(const GObservations& obs) :
          ctool(CTSELECT_NAME, CTSELECT_VERSION)
{
    // Initialise members
    init_members();

    // Set observations
    m_obs = obs;

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
          ctool(CTSELECT_NAME, CTSELECT_VERSION, argc, argv)
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
ctselect::ctselect(const ctselect& app) : ctool(app)
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
        this->ctool::operator=(app);

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
    this->ctool::free_members();
    this->GApplication::free_members();

    // Initialise members
    this->GApplication::init_members();
    this->ctool::init_members();
    init_members();

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

    // Initialise counters
    int n_observations = 0;

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Write header for observation
        if (logTerse()) {
            std::string header = m_obs[i]->instrument() + " observation";
            if (m_obs[i]->name().length() > 1) {
                header += " \"" + m_obs[i]->name() + "\"";
            }
            if (m_obs[i]->id().length() > 1) {
                header += " (id=" + m_obs[i]->id() +")";
            }
            log.header3(header);
        }

        // Initialise event input and output filenames
        m_infiles.push_back("");

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Skip observation if it's not CTA
        if (obs == NULL) {
            if (logTerse()) {
                log << " Skipping ";
                log << m_obs[i]->instrument();
                log << " observation" << std::endl;
            }
            continue;
        }

        // Skip observation if we have a binned observation
        if (obs->eventtype() == "CountsCube") {
            if (logTerse()) {
                log << " Skipping binned ";
                log << obs->instrument();
                log << " observation" << std::endl;
            }
            continue;
        }

        // Increment counter
        n_observations++;

        // Save event file name (for possible saving)
        m_infiles[i] = obs->eventfile();

        // Fall through in case that the event file is empty
        if (obs->events()->size() == 0) {
            if (logTerse()) {
                log << " Warning: No events in event file \"";
                log << m_infiles[i] << "\". Event selection skipped.";
                log << std::endl;
            }
            continue;
        }

        // Get temporary file name
        #if G_USE_MKSTEMP
        char tpl[]  = "ctselectXXXXXX";
        int  fileid = mkstemp(tpl);
        std::string filename(tpl);
        #else
        std::string filename = std::tmpnam(NULL);
        #endif

        // Save observation in temporary file
        obs->save(filename, true);

        // Log saved FITS file
        if (logExplicit()) {
            GFits tmpfile(filename);
            log << std::endl;
            log.header1("FITS file content of temporary file");
            log << tmpfile << std::endl;
            tmpfile.close();
        }

        // Check temporary file
        std::string message = check_infile(filename);
        if (message.length() > 0) {
            throw GException::app_error(G_RUN, message);
        }

        // Load observation from temporary file, including event selection
        select_events(obs, filename);

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
    if (logTerse()) {
        log << std::endl;
        log.header1("Observations after selection");
        log << m_obs << std::endl;
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
    // Write header
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Save observations");
        }
        else {
            log.header1("Save observation");
        }
    }

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
    m_ra     = -1.0;
    m_dec    = -1.0;
    m_rad    = -1.0;
    m_tmin   = 0.0;
    m_tmax   = 0.0;
    m_emin   = 0.0;
    m_emax   = 0.0;
    m_expr.clear();
    m_usethres.clear();

    // Initialise protected members
    m_obs.clear();
    m_infiles.clear();
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
    m_outobs  = app.m_outobs;
    m_prefix   = app.m_prefix;
    m_usepnt   = app.m_usepnt;
    m_ra       = app.m_ra;
    m_dec      = app.m_dec;
    m_rad      = app.m_rad;
    m_tmin     = app.m_tmin;
    m_tmax     = app.m_tmax;
    m_emin     = app.m_emin;
    m_emax     = app.m_emax;
    m_expr     = app.m_expr;
    m_usethres = app.m_usethres;

    // Copy protected members
    m_obs           = app.m_obs;
    m_infiles       = app.m_infiles;
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

    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {

        // Throw exception if no input observation file is given
        require_inobs(G_GET_PARAMETERS);

        // Get observation container without response (not needed)
        m_obs = get_observations(false);

    } // endif: there was no observation in the container

    // Get parameters
    m_usepnt = (*this)["usepnt"].boolean();
    if (!m_usepnt) {

        // Check RA/DEC parameters for validity to read
        if ((*this)["ra"].is_valid() && (*this)["dec"].is_valid()) {
            m_ra         = (*this)["ra"].real();
            m_dec        = (*this)["dec"].real();
            m_select_roi = true;
        }
        else {
            m_select_roi = false;
        }
    }

    // Check if radius is vaild for a RoI selection
    if (m_select_roi && (*this)["rad"].is_valid()) {
        m_rad = (*this)["rad"].real();
    }
    else {
        m_select_roi = false;
    }

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

    // Check parameters (Kluge as the parfile interface does not yet check
    // parameters; this can be removed once this is done properly)
    if (m_usethres != "NONE" && m_usethres != "DEFAULT" && m_usethres != "USER") {
        std::string msg = "Parameter \"usethres\" must be one of \"NONE\", "
                          "\"DEFAULT\" or \"USER\".";
        throw GException::invalid_value(G_GET_PARAMETERS, msg);
    }

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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Select events
 *
 * @param[in] obs CTA observation.
 * @param[in] filename File name.
 *
 * @exception GException::invalid_value
 *            No EVENTS extension found in FITS file.
 *
 * Select events from a FITS file by making use of the selection possibility
 * of the cfitsio library on loading a file. A selection string is created
 * from the specified criteria that is appended to the filename so that
 * cfitsio will automatically filter the event data. This selection string
 * is then applied when opening the FITS file. The opened FITS file is then
 * saved into a temporary file which is the loaded into the actual CTA
 * observation, overwriting the old CTA observation. The ROI, GTI and Ebounds
 * of the CTA event list are then set accordingly to the specified selection.
 * Finally, the temporary file created during this process is removed.
 *
 * Good Time Intervals of the observation will be limited to the time
 * interval [m_tmin, m_tmax]. If m_tmin=m_tmax=0, no time selection is
 * performed.
 ***************************************************************************/
void ctselect::select_events(GCTAObservation* obs, const std::string& filename)
{
    // Initialise selection string
    std::string selection;

    // Initialise selection flags
    bool remove_all = false;

    // Get CTA event list pointer
    GCTAEventList* list =
        static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));

    // Get existing Roi and energy bounds for possible later use
    // (will be empty if unavailable)
    GCTARoi old_roi = list->roi();
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
    bool select_expr   = (gammalib::strip_whitespace(m_expr).length() > 0);

    // Set RA/DEC selection. If the "usepnt" parameter is set to true then
    // use the pointing direction as the ROI centre.
    double ra  = m_ra;
    double dec = m_dec;
    double rad = m_rad;
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
        if (logTerse()) {
            log << gammalib::parformat("Time range");
            log << tmin << " - " << tmax << " s" << std::endl;
        }

    } // endif: made time selection

    // Make energy selection
    if (m_select_energy) {

        // If we have aready a selection than add an "&&" operator
        if (selection.length() > 0) {
            selection += " && ";
        }

        // Log the requested energy selection
        if (logTerse()) {
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
                selection += "ENERGY >= "+std::string(cmin)+" && ENERGY <= "+
                             std::string(cmax);
            }
        }
        else if (emin > 0.0) {
            char cmin[80];
            sprintf(cmin, "%.8f", emin);
            selection += "ENERGY >= "+std::string(cmin);
        }
        else if (emax > 0.0) {
            char cmax[80];
            sprintf(cmax, "%.8f", emax);
            selection += "ENERGY <= "+std::string(cmax);
        }
        else {
            log << "None" << std::endl;
            remove_all = true;
        }

    } // endif: made energy selection

    // Make ROI selection
    if (m_select_roi) {

        // If we have aready a selection than add an "&&" operator
        if (selection.length() > 0) {
            selection += " && ";
        }

        // Log the requested ROI
        if (logTerse()) {
            log << gammalib::parformat("Requested ROI");
            log << "Centre(RA,DEC)=(" << ra << ", " << dec << ") deg, ";
            log << "Radius=" << m_rad << " deg" << std::endl;
        }

        // If we have already an ROI then make sure that the selected
        // ROI overlaps with the existing ROI
        double roi_radius = list->roi().radius();
        if (roi_radius > 0.0) {
            GSkyDir roi_centre = list->roi().centre().dir();
            if (logTerse()) {
                log << gammalib::parformat("ROI of data");
                log << "Centre(RA,DEC)=(" << roi_centre.ra_deg() << ", ";
                log << roi_centre.dec_deg() << ") deg, ";
                log << "Radius=" << roi_radius << " deg" << std::endl;
            }
            GSkyDir centre;
            centre.radec_deg(ra, dec);
            double distance = centre.dist_deg(roi_centre);
            if (distance + rad > roi_radius) {
                rad = roi_radius - distance;
            }
        }

        // Log the selected ROI
        if (logTerse()) {
            log << gammalib::parformat("Selected ROI");
            if (rad <= 0.0) {
                log << "None. There is no overlap between existing ";
                log << "and requested ROI." << std::endl;
            }
            else {
                log << "Centre(RA,DEC)=(" << ra << ", " << dec << ") deg, ";
                log << "Radius=" << rad << " deg" << std::endl;
            }
        }
        
        // Format ROI selection
        char cra[80];
        char cdec[80];
        char crad[80];
        sprintf(cra,  "%.6f", ra);
        sprintf(cdec, "%.6f", dec);
        sprintf(crad, "%.6f", rad);
        selection += "ANGSEP("+std::string(cra)+"," +
                     std::string(cdec)+",RA,DEC) <= " +
                     std::string(crad);

    } // endif: made ROI selection

    // Make an expression selection
    if (select_expr) {

        // If we have aready a selection than add an "&&" operator
        if (selection.length() > 0) {
            selection += " && ";
        }

        // Append the expression
        selection += "("+gammalib::strip_whitespace(m_expr)+")";

    } // endif: made expression selection

    // Dump cfitsio selection string
    if (logTerse()) {
        log << gammalib::parformat("cfitsio selection");
        log << selection << std::endl;
    }

    // Build input filename including selection expression
    std::string expression = filename;
    if (selection.length() > 0) {
        expression += "[EVENTS]["+selection+"]";
    }

    // Dump FITS filename including selection expression
    if (logTerse()) {
        log << gammalib::parformat("FITS filename");
        log << expression << std::endl;
    }

    // Open FITS file
    GFits file(expression);

    // Log selected FITS file
    if (logExplicit()) {
        log << std::endl;
        log.header1("FITS file content after selection");
        log << file << std::endl;
    }

    // Check if we have an EVENTS HDU
    if (!file.contains("EVENTS")) {
        std::string msg = "No \"EVENTS\" extension found in FITS file \""+
                          expression+"\".";
        throw GException::invalid_value(G_SELECT_EVENTS, msg);
    }

    // Determine number of events in EVENTS HDU
    int nevents = file.table("EVENTS")->nrows();

    // If the selected event list is empty or if removal of all events
    // has been requested then append an empty event list to the observation.
    if ((nevents < 1) || (remove_all)) {

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

    // If ROI selection has been applied then set the event list ROI
    if (m_select_roi) {
        GCTAInstDir instdir;
        instdir.dir().radec_deg(ra, dec);
        list->roi(GCTARoi(instdir, rad));
    } // endif: Roi selection was performed

    else if (old_roi.is_valid()) {
        // Restore old Roi information in case no selection was performed
        // and RoI was existing before
        list->roi(old_roi);
    }

    // Set event list GTI (in any case as any event list has a GTI)
    list->gti(gti);

    // If an energy selection has been applied then set the energy boundaries
    if (m_select_energy) {
        GEbounds ebounds;
        ebounds.append(GEnergy(emin, "TeV"), GEnergy(emax, "TeV"));
        list->ebounds(ebounds);
    } //endif: energy selection was performed

    else if (old_ebounds.size() > 0) {
        // Restore old Ebounds in case no energy selection was performed
        // and observation already had valid Ebounds
        list->ebounds(old_ebounds);
    }

    // Recompute ontime and livetime.
    GTime meantime = 0.5 * (list->gti().tstart() + list->gti().tstop());
    obs->ontime(list->gti().ontime());
    obs->livetime(list->gti().ontime() * obs->deadc(meantime));

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
        const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(obs->response());
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
        if (missing.size() > 0) {
            message = "The following columns are missing in the"
                      " \"EVENTS\" extension of input file \""
                    + m_outobs + "\": ";
            for (int i = 0; i < missing.size(); ++i) {
                message += "\"" + missing[i] + "\"";
                if (i < missing.size()-1) {
                    message += ", ";
                }
            }
        }

    }
    catch (GException::fits_hdu_not_found& e) {
        message = "No \"EVENTS\" extension found in input file \""
                + m_outobs + "\".";
    }

    // Return
    return message;
}


/***********************************************************************//**
 * @brief Set output file name.
 *
 * @param[in] filename Input file name.
 *
 * Converts an input filename into an output filename by prepending the
 * prefix stored in the member m_prefix to the input filename. Any path will
 * be stripped from the input filename. Also a trailing ".gz" will be
 * stripped.
 ***************************************************************************/
std::string ctselect::set_outfile_name(const std::string& filename) const
{
    // Split input filename into path elements
    std::vector<std::string> elements = gammalib::split(filename, "/");

    // The last path element is the filename
    std::string outname = m_prefix + elements[elements.size()-1];

    // Strip any ".gz"
    outname = gammalib::strip_chars(outname, ".gz");
    
    // Return output filename
    return outname;
}


/***********************************************************************//**
 * @brief Save event list in FITS format.
 *
 * Save the event list as a FITS file. The filename of the FITS file is
 * specified by the outfile parameter.
 ***************************************************************************/
void ctselect::save_fits(void)
{
    // Get output filename
    m_outobs = (*this)["outobs"].filename();

    // Get CTA observation from observation container
    GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

    // Save only if it's a CTA observation
    if (obs != NULL) {
    
        // Save only if filename is non-empty
        if (m_infiles[0].length() > 0) {

            // Dump filename
            if (logTerse()) {
                log << "Save \""+m_infiles[0]+"\"" << std::endl;
            }

            // Save event list
            save_event_list(obs, m_infiles[0], m_outobs);

        } // endif: filename was non empty

    } // endif: observation was CTA observation

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
    m_prefix  = (*this)["prefix"].string();

    // Issue warning if output filename has no .xml suffix
    std::string suffix = gammalib::tolower(m_outobs.substr(m_outobs.length()-4,4));
    if (suffix != ".xml") {
        log << "*** WARNING: Name of observation definition output file \""+
               m_outobs+"\"" << std::endl;
        log << "*** WARNING: does not terminate with \".xml\"." << std::endl;
        log << "*** WARNING: This is not an error, but might be misleading."
               " It is recommended" << std::endl;
        log << "*** WARNING: to use the suffix \".xml\" for observation"
               " definition files." << std::endl;
    }

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

        // Dump filename
        if (logTerse()) {
            log << "Save \""+outfile+"\"" << std::endl;
        }

        // Store output file name in observation
        obs->eventfile(outfile);

        // Save event list
        save_event_list(obs, m_infiles[i], outfile);

    } // endfor: looped over observations

    // Dump filename
    if (logTerse()) {
        log << "Save \""+m_outobs+"\"" << std::endl;
    }

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
 * @param[in] outfile Output file name.
 *
 * Saves an event list into a FITS file and copy all others extensions from
 * the input file to the output file.
 ***************************************************************************/
void ctselect::save_event_list(const GCTAObservation* obs,
                               const std::string&     infile,
                               const std::string&     outfile) const
{
    // Save only if observation is valid
    if (obs != NULL) {

        // Save only if we have an event list
        if (obs->eventtype() == "EventList") {

            // Save observation into FITS file
            obs->save(outfile, clobber());

            // Copy all extensions other than EVENTS and GTI from the input to
            // the output event list. The EVENTS and GTI extensions are written
            // by the save method, all others that may eventually be present
            // have to be copied by hand.
            GFits infits(infile);
            GFits outfits(outfile);
            for (int extno = 1; extno < infits.size(); ++extno) {
                GFitsHDU* hdu = infits.at(extno);
                if (hdu->extname() != "EVENTS" && hdu->extname() != "GTI") {
                    outfits.append(*hdu);
                }
            }

            // Close input file
            infits.close();

            // Save file to disk and close it (we need both operations)
            outfits.save(true);
            outfits.close();

        } // endif: observation was unbinned

    } // endif: observation was valid

    // Return
    return;
}
