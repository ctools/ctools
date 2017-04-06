/***************************************************************************
 *           ctprob - Computes probability for a given model               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Leonardo Di Venere                          *
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
 * @file ctprob.cpp
 * @brief Computes probability for a given model
 * @author Leonardo Di Venere
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>
#include "ctprob.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_RUN                                               "ctprob::run()"
#define G_GET_OBS                                       "ctprob::get_obs()"

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
ctprob::ctprob(void) : ctobservation(CTPROB_NAME, CTPROB_VERSION)
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
ctprob::ctprob(const GObservations& obs) :
          ctobservation(CTPROB_NAME, CTPROB_VERSION, obs)
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
ctprob::ctprob(int argc, char *argv[]) : 
          ctobservation(CTPROB_NAME, CTPROB_VERSION, argc, argv)
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
ctprob::ctprob(const ctprob& app) : ctobservation(app)
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
ctprob::~ctprob(void)
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
ctprob& ctprob::operator=(const ctprob& app)
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
 * @brief Clear ctprob tool
 *
 * Clears ctprob tool.
 ***************************************************************************/
void ctprob::clear(void)
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
 * @brief Append probability columns event data
 *
 * This method reads in the application parameters and loops over all
 * observations that were found to perform the probability calculation.
 * Each observation is written to a temporary file, which is then 
 * re-opened and used for the actual calculation. The temporary file
 *  is deleted after this action so that no disk overflow will occur.
 ***************************************************************************/
void ctprob::run(void)
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
        char tpl[]  = "ctprobXXXXXX";
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
        evaluate_probability(obs);
	//select_events(obs, filename, m_evtname[i], m_gtiname[i]);

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
void ctprob::save(void)
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
void ctprob::publish(const std::string& name)
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
                    user_name = CTPROB_NAME;
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
void ctprob::init_members(void)
{
    // Initialise parameters
    m_outobs.clear();
    m_prefix.clear();
    /*    m_usepnt = false;
    m_roi.clear();
    m_tmin   = 0.0;
    m_tmax   = 0.0;
    m_emin   = 0.0;
    m_emax   = 0.0;
    m_expr.clear();
    m_usethres.clear();*/
    m_chatter = static_cast<GChatter>(2);

    // Initialise protected members
    m_infiles.clear();
    m_evtname.clear();
    m_gtiname.clear();
    //    m_timemin.clear();
    //m_timemax.clear();
    //m_select_energy = true;
    //m_select_roi    = true;
    //m_select_time   = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctprob::copy_members(const ctprob& app)
{
    // Copy parameters
    m_outobs   = app.m_outobs;
    m_prefix   = app.m_prefix;
    /*    m_usepnt   = app.m_usepnt;
    m_roi      = app.m_roi;
    m_tmin     = app.m_tmin;
    m_tmax     = app.m_tmax;
    m_emin     = app.m_emin;
    m_emax     = app.m_emax;
    m_expr     = app.m_expr;
    m_usethres = app.m_usethres;*/
    m_chatter  = app.m_chatter;

    // Copy protected members
    m_infiles       = app.m_infiles;
    m_evtname       = app.m_evtname;
    m_gtiname       = app.m_gtiname;
    /*    m_timemin       = app.m_timemin;
    m_timemax       = app.m_timemax;
    m_select_energy = app.m_select_energy;
    m_select_roi    = app.m_select_roi;
    m_select_time   = app.m_select_time;*/

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctprob::free_members(void)
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
void ctprob::get_parameters(void)
{
    // Initialise selection flags
  /*    m_select_energy = true;
    m_select_roi    = true;
    m_select_time   = true;
  */
    // Setup observations from "inobs" parameter. Do not request response
    // information and do not accept counts cubes.
    if (m_obs.size() == 0) {
        get_obs();
    }
    // ... otherwise add response information and energy boundaries in case
    // that they are missing
    else {
        setup_observations(m_obs);
    }
    //setup_observations(m_obs, false, true, false);

    // Read model definition file if required
    if (m_obs.models().size() == 0) {

        // Get model filename
        std::string inmodel = (*this)["inmodel"].filename();

        // Load models from file
        m_obs.models(inmodel);

    } // endif: there were no models

    // Get energy dispersion flag parameters
    m_apply_edisp = (*this)["edisp"].boolean();


    // Get remaining parameters
    m_publish = (*this)["publish"].boolean();
    m_chatter = static_cast<GChatter>((*this)["chatter"].integer());

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
 * @brief Get observation container
 *
 * @exception GException::invalid_value
 *            Invalid FITS file.
 *
 * Get an observation container according to the user parameters. The method
 * supports loading of a individual FITS file or an observation definition
 * file in XML format.
 ***************************************************************************/
void ctprob::get_obs(void)
{
    // Get the filename from the input parameters
    std::string filename = (*this)["inobs"].filename();

    // If no observation definition file throw an exception 
    if (!is_valid_filename(filename)) {
      throw GException::invalid_value(G_GET_OBS, "Input event list is not valid.");
    } // endif: filename was not valid

    // ... otherwise we have a file name
    else {

        // If file is a FITS file then create an empty CTA observation
        // and load file into observation
        if (GFilename(filename).is_fits()) {

            // Allocate empty CTA observation
            GCTAObservation cta;

            // Load data
            cta.load(filename);

            // Set response
            set_obs_response(&cta);

            // Append observation to container
            m_obs.append(cta);

            // Signal that no XML file should be used for storage
            m_use_xml = false;

        }

        // ... otherwise load file into observation container
        else {

            // Load observations from XML file
            m_obs.load(filename);

            // For all observations that have no response, set the response
            // from the task parameters
            set_response(m_obs);

            // Set observation boundary parameters (emin, emax, rad)
            set_obs_bounds(m_obs);

            // Signal that XML file should be used for storage
            m_use_xml = true;

        } // endelse: file was an XML file

    }

    // Return
    return;

}



/***********************************************************************//**
 * @brief Evaluate probability for events
 *
 * @param[in,out] obs CTA observation.
 *
 *
 * For each event from the input event file evaluates the probability that
 * the event comes from any of the source in the input model. This is done 
 * by evaluating differential expected counts for the event direction and 
 * energy and normalizing it to the total differential expected counts for 
 * given model.
 * observation is replaced by selected event list read from the FITS file.
 *
 * A FITS column for each source in the model is created and added to 
 * the event list. The name of each column is made by the source name with 
 * the prefix "PROB_".
 ***************************************************************************/
void ctprob::evaluate_probability(GCTAObservation*   obs)
{
    // Write header into logger
    log_header3(NORMAL, "Evaluating probability for events");

    GCTAEventList* evts = static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));
    const GEvent* evt;

    double total=0.;
    double value;

    std::vector<GFitsTableFloatCol*> columns;

    //define colums
    for (int j = 0; j<m_obs.models().size(); ++j) {
      const std::string mdl_name = m_obs.models()[j] ->name();
      GFitsTableFloatCol* col = new GFitsTableFloatCol("PROB_"+mdl_name, evts->size());
      columns.push_back(col);
    }

    // Loop over events
    for (int i = 0; i < evts->size(); ++i) {
      total = 0.;
      evt = (*evts)[i];
      std::vector<double> values;
      //Loop over models
      for (int j = 0; j<m_obs.models().size(); ++j) {
	value =  m_obs.models()[j]->eval( *evt, *obs);
	values.push_back(value);
	total += value;
      }
      for (int j = 0; j<m_obs.models().size(); ++j) {
	values[j] /= total;
	(*(columns[j]))(i) = values[j];
      }
    }

    //append columns to observation event list
    for (int j = 0; j<m_obs.models().size(); ++j) {
      evts->append_column(*columns[j]);
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
std::string ctprob::check_infile(const std::string& filename,
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
std::string ctprob::set_outfile_name(const std::string& filename) const
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
std::string ctprob::get_gtiname(const std::string& filename,
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
void ctprob::save_fits(void)
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
void ctprob::save_xml(void)
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
void ctprob::save_event_list(const GCTAObservation* obs,
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
