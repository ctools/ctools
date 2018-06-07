/***************************************************************************
 *          ctprob - Computes event probability for a given model          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2018 by Leonardo Di Venere                          *
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
 * @brief Event probability computation tool implementation
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
#define G_EVALUATE_PROB      "ctprob::evaluate_probability(GCTAObservation*)"

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
ctprob::ctprob(void) : ctobservation(CTPROB_NAME, VERSION)
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
        ctobservation(CTPROB_NAME, VERSION, obs)
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
        ctobservation(CTPROB_NAME, VERSION, argc, argv)
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
 * @brief Runs the ctprob tool
 *
 * This method reads in the application parameters and loops over all
 * unbinned observations to compute the probability for each event that it
 * comes from one of the model components.
 ***************************************************************************/
void ctprob::run(void)
{
    // Switch screen logging on in debug mode
    if (logDebug()) {
        log.cout(true);
    }

    // Get parameters
    get_parameters();

    // Set energy dispersion flags of all CTA observations and save old
    // values in save_edisp vector
    std::vector<bool> save_edisp = set_edisp(m_obs, m_apply_edisp);

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write header into logger
    log_header1(TERSE, "Compute event probabilities");

    // Initialise counters
    int n_observations = 0;

    // Loop over all unbinned CTA observations in the container
    for (GCTAObservation* obs = first_unbinned_observation(); obs != NULL;
         obs = next_unbinned_observation()) {

        // Increment counter
        n_observations++;

        // Evaluate event probabilities
        evaluate_probability(obs);

    } // endfor: looped over observations

    // If more than a single observation has been handled then make sure that
    // an XML file will be used for storage
    if (n_observations > 1) {
        m_use_xml = true;
    }

    // Write observation(s) into logger
    log_observations(NORMAL, m_obs, "Output observation");

    // Restore energy dispersion flags of all CTA observations
    restore_edisp(m_obs, save_edisp);

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
    m_apply_edisp = false;
    m_publish     = false;
    m_chatter     = static_cast<GChatter>(2);

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
    m_apply_edisp = app.m_apply_edisp;
    m_publish     = app.m_publish;
    m_chatter     = app.m_chatter;

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
 * Get all application parameters by either querying the parameters or by
 * retrieving the parameters for the parameter file.
 ***************************************************************************/
void ctprob::get_parameters(void)
{
    // Setup observations from "inobs" parameter. Request response
    // information, accept event lists but do not accept counts cubes.
    setup_observations(m_obs, true, true, false);

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

    // If needed later, query output filename and prefix now
    if (read_ahead()) {
        (*this)["outobs"].query();
        (*this)["prefix"].query();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Evaluate probability for events
 *
 * @param[in,out] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            Observation does not contain an event list.
 *
 * Computes for all events the probability that an event arises from a
 * specific model component. The computation is done by by evaluating the
 * differential event probability for a given model and normalizing these
 * probabilities so that the sum of the probabilities for all model
 * components is unity for each event.
 *
 * The method appends for each model component a single precision column to
 * the event list that contains the event probabilities. The name of the
 * column is build from the model name prefixed with `PROB_`.
 ***************************************************************************/
void ctprob::evaluate_probability(GCTAObservation* obs)
{
    // Write header into logger
    log_header3(NORMAL, "Evaluate events probabilities");

    // Make sure that the observation holds a CTA event list. If this
    // is not the case then throw an exception.
    GCTAEventList* events = dynamic_cast<GCTAEventList*>(const_cast<GEvents*>
                            (obs->events()));
    if (events == NULL) {
        std::string msg = "CTA Observation does not contain an event list. An "
                          "event list is needed to compute the event "
                          "probabilities.";
        throw GException::invalid_value(G_EVALUATE_PROB, msg);
    }

    // Continue only if there are model components
    int nmodels = m_obs.models().size();
    if (nmodels > 0) {

        // Allocate single precision FITS table columns for all model components
        std::vector<GFitsTableFloatCol*> columns;
        for (int k = 0; k < nmodels; ++k) {
            const std::string   colname = "PROB_" + m_obs.models()[k] ->name();
            GFitsTableFloatCol* col     =
                                new GFitsTableFloatCol(colname, events->size());
            columns.push_back(col);
        }

        // Loop over all events in the event list
        for (int i = 0; i < events->size(); ++i) {

            // Initialise total probability
            double total = 0.0;

            // Get point to event
            const GEvent* event = (*events)[i];

            // Initialise probability vector
            std::vector<double> probabilities(nmodels, 0.0);

            // Loop over models
            for (int k = 0; k < nmodels; ++k) {
                double value      =  m_obs.models()[k]->eval(*event, *obs);
                probabilities[k]  = value;
                total            += value;
            }

            // Normalize probabilities to unity
            if (total > 0.0) {
                for (int k = 0; k < nmodels; ++k) {
                    probabilities[k] /= total;
                }
            }

            // Store probabilities in FITS columns
            for (int k = 0; k < nmodels; ++k) {
                (*(columns[k]))(i) = probabilities[k];
            }
	
        } // endfor: looped over all events

        // Append columns to event list
        for (int k = 0; k < nmodels; ++k) {
            events->append_column(*columns[k]);
        }

    } // endif: there were model components

    // Return
    return;
}
