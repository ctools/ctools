/***************************************************************************
 *          ctphase - Append phase information to CTA events file          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Joshua Cardenzana                                *
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
 * @author Joshua Cardenzana
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>
#include "ctphase.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                          "ctphase::get_parameters()"
#define G_PHASE_EVENTS              "ctphase::phase_events(GCTAObservation*)"

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
ctphase::ctphase(void) : ctobservation(CTPHASE_NAME, VERSION)
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
         ctobservation(CTPHASE_NAME, VERSION, obs)
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
         ctobservation(CTPHASE_NAME, VERSION, argc, argv)
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
 * @brief Runs the ctprob tool
 *
 * This method reads in the application parameters and loops over all
 * unbinned observations to compute the phase information for each event
 * based on an input model.
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
    log_header1(TERSE, "Compute event phase");

    // Initialise counters
    int n_observations = 0;

    // Loop over all unbinned CTA observations in the container
    for (GCTAObservation* obs = first_unbinned_observation(); obs != NULL;
         obs = next_unbinned_observation()) {

        // Increment counter
        n_observations++;

        // Compute event phase
        phase_events(obs);

    } // endfor: looped over observations

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
        save_events_xml();
    }

    // Case B: Save event file as FITS file
    else {
        save_events_fits();
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
    m_chatter = static_cast<GChatter>(2);

    // Initialise protected members
    m_phase.clear();
    
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
    m_chatter = app.m_chatter;

    // Copy protected members
    m_phase = app.m_phase;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctphase::free_members(void)
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
void ctphase::get_parameters(void)
{

    // Setup observations from "inobs" parameter. Do not request response
    // information and do not accept counts cubes.
    setup_observations(m_obs, false, true, false);

    // If a model definition XML file is provided the
    GFilename inmodel = (*this)["inmodel"].filename();
    if (is_valid_filename(inmodel)) {

        // Load the models from the XML file
        GModels models(inmodel);

        // Get the source component
        std::string source = (*this)["source"].string();

        // Throw an exception if source is not found in the model container
        if (!models.contains(source)) {
            std::string msg = "Source model \""+source+"\" not found in model "
                              "definition XML file \""+inmodel.url()+"\".";
            throw GException::invalid_value(G_GET_PARAMETERS, msg);
        }

        // Get pointer to source model
        GModelSky* sky = dynamic_cast<GModelSky*>(models[source]);

        // Throw an exception if model is not a sky model
        if (sky == NULL) {
            std::string msg = "Source model \""+source+"\" is not a sky model.";
            throw GException::invalid_value(G_GET_PARAMETERS, msg);
        }

        // Get pointer to temporal phase curve
        GModelTemporalPhaseCurve* phasecurve =
              dynamic_cast<GModelTemporalPhaseCurve*>(sky->temporal());

        // Throw an exception if model has not a temporal phase curve
        // component
        if (phasecurve == NULL) {
            std::string msg = "Source model \""+source+"\" does not have a "
                              "temporal component with phase information.";
            throw GException::invalid_value(G_GET_PARAMETERS, msg);
        }

        // Set phase curve
        m_phase = GModelTemporalPhaseCurve(*phasecurve);

    } // endif: model definition XML file provided

    // Otherwise read phase curve information from the parameter files
    else {

        // Read phase curve user parameters
        double mjd_value = (*this)["mjd"].real();
        double phase     = (*this)["phase"].real();
        double f0        = (*this)["f0"].real();
        double f1        = (*this)["f1"].real();
        double f2        = (*this)["f2"].real();

        // Set reference time
        GTime mjd;
        mjd.mjd(mjd_value);

        // Set phase curve
        m_phase = GModelTemporalPhaseCurve();
        m_phase.mjd(mjd);
        m_phase.phase(phase);
        m_phase.f0(f0);
        m_phase.f1(f1);
        m_phase.f2(f2);

    } // endelse: read phase curve information from User parameters

    // Get other User parameters
    m_chatter = static_cast<GChatter>((*this)["chatter"].integer());

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        (*this)["outobs"].filename();
        (*this)["prefix"].string();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute event phase for an observation
 *
 * @param[in,out] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No events extension found in FITS file.
 *
 * Append phase information to events from a FITS file by adding a "PHASE"
 * column to the FITS file.
 ***************************************************************************/
void ctphase::phase_events(GCTAObservation* obs)
{
    // Write header into logger
    log_header3(NORMAL, "Compute event phases");

    // Make sure that the observation holds a CTA event list. If this
    // is not the case then throw an exception.
    GCTAEventList* events = dynamic_cast<GCTAEventList*>(const_cast<GEvents*>
                            (obs->events()));
    if (events == NULL) {
        std::string msg = "CTA Observation does not contain an event list. An "
                          "event list is needed to compute the event phases.";
        throw GException::invalid_value(G_PHASE_EVENTS, msg);
    }

    // Continue only if there are events
    if (events->size() > 0) {

        // Print a warning into the log file if there are some precision concerns
        GCTAEventAtom* event = (*events)[0];
        double         dt    = std::abs(event->time() - m_phase.mjd());
        if (((m_phase.f1() * dt) > 1.0e8) ||
             (m_phase.f2() * dt * dt) > 1.0e8) {
            log_string(NORMAL,"*** WARNING ***");
            log_string(NORMAL,"   Values supplied for reference MJD and/or f1 and/or f2");
            log_string(NORMAL,"   are large and may result in numerical precision issues.");
            log_string(NORMAL,"   Consider using an ephemeris derived from a date closer");
            log_string(NORMAL,"   to the data being analyzed.");
            log_value(NORMAL, "     Observation ID", obs->obs_id());
            log_value(NORMAL, "     Reference MJD", m_phase.mjd().mjd());
            log_value(NORMAL, "     First event MJD", event->time().mjd());
            log_string(NORMAL,"*** WARNING ***");
        }

        // Compute the phase for all events
        for (int i = 0; i < events->size(); ++i) {

            // Get pointer to the event
            GCTAEventAtom* event = (*events)[i];
        
            // Get the event time
            GTime time = event->time();
        
            // Set the event phase
            event->phase(m_phase.phase(time));

        } // endfor: looped over events

    } // endif: there were events
    
    // Signal that phase information is available
    events->has_phase(true);

    // Return
    return;
}
