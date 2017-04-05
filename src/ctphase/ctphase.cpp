/***************************************************************************
 *          ctphase - Append phase information to CTA events file          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Leonardo Di Venere                               *
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
 * @file ctphase.cpp
 * @brief Append phase information to CTA events file
 * @author Leonardo Di Venere
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>
#include "ctphase.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_RUN                                               "ctphase::run()"
#define G_PHASE_EVENTS          "ctphase::phase_events(GCTAObservation*, "\
                                  "std::string&, std::string&, std::string&)"
#define G_SET_EBOUNDS    "ctphase::set_ebounds(GCTAObservation*, GEbounds&)"

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
ctphase::ctphase(void) : ctobservation(CTPHASE_NAME, CTPHASE_VERSION)
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
ctphase::ctphase(const GObservations& obs) :
          ctobservation(CTPHASE_NAME, CTPHASE_VERSION, obs)
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
ctphase::ctphase(int argc, char *argv[]) : 
          ctobservation(CTPHASE_NAME, CTPHASE_VERSION, argc, argv)
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
ctphase::ctphase(const ctphase& app) : ctobservation(app)
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
ctphase::~ctphase(void)
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
ctphase& ctphase::operator=(const ctphase& app)
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
 * @brief Clear ctphase tool
 *
 * Clears ctphase tool.
 ***************************************************************************/
void ctphase::clear(void)
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
 * @brief Append phase information to event data
 *
 * This method reads in the application parameters and loops over all
 * observations that were found to append phase information to the data.
 ***************************************************************************/
void ctphase::run(void)
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
    log_header1(TERSE, "Appending phase informaiton");

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
                       m_infiles[i]+"\". Event phasing skipped.");
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
        char tpl[]  = "ctphaseXXXXXX";
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

        // Load observation from temporary file and compute phase
        phase_events(obs, filename, m_evtname[i], m_gtiname[i]);

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
 * @brief Save the output event list(s)
 *
 * This method saves the updated event list(s) into FITS file(s). There are
 * two modes, depending on the m_use_xml flag.
 *
 * If m_use_xml is true, all output event list(s) will be saved into FITS
 * files, where the output filenames are constructued from the input
 * filenames by prepending the m_prefix string to name. Any path information
 * will be stripped form the input name, hence event files will be written
 * into the local working directory (unless some path information is present
 * in the prefix). In addition, an XML file will be created that gathers
 * the filename information for the output event list(s). If an XML file
 * was present on input, all metadata information will be copied from this
 * input file.
 *
 * If m_use_xml is false, the output event list will be saved into a FITS
 * file.
 ***************************************************************************/
void ctphase::save(void)
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
void ctphase::publish(const std::string& name)
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
                    user_name = CTPHASE_NAME;
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
void ctphase::init_members(void)
{
    // Initialise parameters
    m_outobs.clear();
    m_prefix.clear();
    m_phase.clear();
    m_chatter = static_cast<GChatter>(2);

    // Initialise protected members
    m_infiles.clear();
    m_evtname.clear();
    m_gtiname.clear();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctphase::copy_members(const ctphase& app)
{
    // Copy parameters
    m_outobs   = app.m_outobs;
    m_prefix   = app.m_prefix;
    m_chatter  = app.m_chatter;

    // Copy protected members
    m_infiles       = app.m_infiles;
    m_evtname       = app.m_evtname;
    m_gtiname       = app.m_gtiname;
    m_phase         = app.m_phase;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctphase::free_members(void)
{
    m_phase.clear();
    
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
void ctphase::get_parameters(void)
{

    // Setup observations from "inobs" parameter. Do not request response
    // information and do not accept counts cubes.
    setup_observations(m_obs, false, true, false);

    // Check for sanity of frequency and phase information
    if ((*this)["p0"].is_valid() &&
        (*this)["mjd"].is_valid() &&
        (*this)["f0"].is_valid()) {
        
        // The following creates an XML object responsible for initializing
        // 'm_phase'. This was done as the default values for the parameters in
        // 'GModelTemporalPhaseCurve' are such that the ranges are predefined.
        // Initializing 'm_phase' from an XML object circumvents this restriction.
        
        GXmlElement norm;
        norm.name("parameter");
        norm.attribute("name", "Normalization");
        norm.attribute("value", "1");
        norm.attribute("min", "0");
        norm.attribute("max", "1000");
        
        GXmlElement mjd;
        mjd.name("parameter");
        mjd.attribute("name", "MJD");
        mjd.attribute("value", (*this)["mjd"].value());
        mjd.attribute("min", (*this)["mjd"].min());
        mjd.attribute("max", (*this)["mjd"].max());
        
        GXmlElement p0;
        p0.name("parameter");
        p0.attribute("name", "Phase");
        p0.attribute("value", (*this)["p0"].value());
        p0.attribute("min", (*this)["p0"].min());
        p0.attribute("max", (*this)["p0"].max());
        
        GXmlElement f0;
        f0.name("parameter");
        f0.attribute("name", "F0");
        f0.attribute("value", (*this)["f0"].value());
        f0.attribute("min", (*this)["f0"].min());
        f0.attribute("max", (*this)["f0"].max());
        
        GXmlElement f1;
        f1.name("parameter");
        f1.attribute("name", "F1");
        f1.attribute("value", (*this)["f1"].value());
        f1.attribute("min", (*this)["f1"].min());
        f1.attribute("max", (*this)["f1"].max());
        
        GXmlElement f2;
        f2.name("parameter");
        f2.attribute("name", "F2");
        f2.attribute("value", (*this)["f2"].value());
        f2.attribute("min", (*this)["f2"].min());
        f2.attribute("max", (*this)["f2"].max());
        
        GXmlElement info;
        info.append(norm);
        info.append(mjd);
        info.append(p0);
        info.append(f0);
        info.append(f1);
        info.append(f2);
        
        info.attribute("file","");
        
        m_phase = GModelTemporalPhaseCurve();
        m_phase.read(info);
    }
    else {
        GException::invalid_value("ctphase::get_parameters()",
                                  "Invalid parameter value for p0, mjd, or f0.");
    }

    // Get other User parameters
    m_chatter  = static_cast<GChatter>((*this)["chatter"].integer());

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outobs = (*this)["outobs"].filename();
        m_prefix = (*this)["prefix"].string();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Define phase for events
 *
 * @param[in,out] obs CTA observation.
 * @param[in] filename File name.
 * @param[in] evtname Event extension name.
 * @param[in] gtiname GTI extension name.
 *
 * @exception GException::invalid_value
 *            No events extension found in FITS file.
 *
 * Append phase information to events from a FITS file by adding a "PHASE"
 * column to the FITS file.
 ***************************************************************************/
void ctphase::phase_events(GCTAObservation*   obs,
                             const std::string& filename,
                             const std::string& evtname,
                             const std::string& gtiname)
{
    // Write header into logger
    log_header3(NORMAL, "Events phasing");

    // Get CTA event list pointer
    GCTAEventList* list =
        static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));
    
    for (int evnt=0; evnt<list->size(); evnt++) {
        // Get the next event
        GCTAEventAtom* event = (*list)[evnt];
        
        // Get the event time
        GTime time = event->time();
        
        // Set the event phase
        event->phase(m_phase.phase(time));
    }
    
    // Make sure the file knows to save the phase information
    if (!list->has_phase()) {
        list->has_phase(true);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check input filename
 *
 * @param[in] filename File name.
 * @param[in] evtname Event extension name.
 *
 * This method checks if the input FITS file is correct.
 ***************************************************************************/
std::string ctphase::check_infile(const std::string& filename,
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
std::string ctphase::set_outfile_name(const std::string& filename) const
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
std::string ctphase::get_gtiname(const std::string& filename,
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

    } // endif: filename was not empty

    // Return GTI name
    return (gtiname);
}


/***********************************************************************//**
 * @brief Save event list in FITS format.
 *
 * Save the event list as a FITS file. The file name of the FITS file is
 * specified by the "outobs" parameter.
 ***************************************************************************/
void ctphase::save_fits(void)
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
void ctphase::save_xml(void)
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
 * will write the updated events into the "EVENTS1" extension and the
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
void ctphase::save_event_list(const GCTAObservation* obs,
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
