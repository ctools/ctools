/***************************************************************************
 *                   ctlike - CTA maximum likelihood tool                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file ctlike.cpp
 * @brief CTA maximum likelihood tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctlike.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                           "ctlike::get_parameters()"

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
ctlike::ctlike(void) : GApplication(CTLIKE_NAME, CTLIKE_VERSION)
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
 * param[in] obs Observation container.
 *
 * This method creates an instance of the class by copying an existing
 * observations container.
 ***************************************************************************/
ctlike::ctlike(const GObservations& obs) :
        GApplication(CTLIKE_NAME, CTLIKE_VERSION)
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
ctlike::ctlike(int argc, char *argv[]) : 
                        GApplication(CTLIKE_NAME, CTLIKE_VERSION, argc, argv)
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
ctlike::ctlike(const ctlike& app) : GApplication(app)
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
ctlike::~ctlike(void)
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
ctlike& ctlike::operator= (const ctlike& app)
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
void ctlike::clear(void)
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
 * This method performs a maximum likelihood analysis of a observation given
 * in an observation container and saves the results in an XML file.
 ***************************************************************************/
void ctlike::execute(void)
{
    // Signal that some parameters should be read ahead
    m_read_ahead = true;

    // Bin the event data
    run();

    // Save the results into XML file
    save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Run maximum likelihood analysis
 *
 * The following analysis steps are performed:
 * 1. Read the parameters (and write them into logger)
 * 2. Load observation
 * 3. Setup models for optimizing
 * 4. Optimize model (and write result into logger)
 ***************************************************************************/
void ctlike::run(void)
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

    // Set energy dispersion flag for all CTA observations and save old
    // values in save_edisp vector
    std::vector<bool> save_edisp;
    save_edisp.assign(m_obs.size(), false);
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            save_edisp[i] = obs->response()->apply_edisp();
            obs->response()->apply_edisp(m_apply_edisp);
        }
    }

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Observations");
        }
        else {
            log.header1("Observation");
        }
        log << m_obs << std::endl;
    }

    // Optimize model parameters using LM optimizer
    optimize_lm();

    // Optionally perform TS computation
    if (m_tscalc) {
    
        // Store original maximum likelihood and models
        double  logL_src    = m_logL;
        GModels models_orig = m_obs.models();
    
        // Store models with free parameters
        std::vector<std::string> free_srcs;
        GModels models = m_obs.models();
        for (int i = 0; i < models.size(); ++i) {
            GModel* model = models[i];
            for (int j = 0; j < model->size(); ++j) {
                GModelPar par = model->at(j);
                if (par.is_free()){
                    free_srcs.push_back(model->name());
                    break;
                }
            }
        }
    
        // Loop over stored models, remove source and refit
        for (int i = 0; i < free_srcs.size(); ++i) {
            models.remove(free_srcs[i]);  
            m_obs.models(models);    
            double logL_nosrc = reoptimize_lm();
            double ts         = 2.0 * (logL_src-logL_nosrc);
            models_orig[free_srcs[i]]->ts(ts);
            models = models_orig;
        }
    
        // Restore best fit values
        m_obs.models(models_orig);

    } // endif: did optional TS computation

    // Compute number of observed events in all observations
    double num_events = 0.0;
    for (int i = 0; i < m_obs.size(); ++i) {
        num_events += m_obs[i]->events()->number();
    }

    // Write results into logger
    if (logTerse()) {
        log << gammalib::parformat("Maximum log likelihood");
        log << gammalib::str(m_logL,3) << std::endl;
        log << gammalib::parformat("Observed events  (Nobs)");
        log << gammalib::str(num_events, 3) << std::endl;
        log << gammalib::parformat("Predicted events (Npred)");
        log << gammalib::str(m_obs.npred(), 3);
        log << " (Nobs - Npred = ";
        log << gammalib::str(num_events-m_obs.npred());
        log << ")" << std::endl;
        log << m_obs.models() << std::endl;
    }

    // Restore energy dispersion flag for all CTA observations
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            obs->response()->apply_edisp(save_edisp[i]);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save results
 *
 * This method saves the fit results in a XML file.
 ***************************************************************************/
void ctlike::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save results");
    }

    // Get output filename
    m_outmdl = (*this)["outmdl"].filename();

    // Write results out as XML model
    if (gammalib::toupper(m_outmdl) != "NONE") {
        m_obs.models().save(m_outmdl);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all required task parameters from the parameter file or (if specified)
 * by querying the user. Observation dependent parameters will only be read
 * if the observation container is actually empty. Observation dependent
 * parameters are:
 * "stat" (statistics to be used for observation),
 * "caldb" (calibration database),
 * "irf" (instrument response function), and
 * "infile" (input file name).
 * The model will only be loaded if no model components exist in the
 * observation container.
 *
 * This method handles both loading of FITS files and of handling XML
 * observation definition files.
 ***************************************************************************/
void ctlike::get_parameters(void)
{
    // If there are no observations in container then add a single CTA
    // observation using the parameters from the parameter file
    if (m_obs.size() == 0) {

        // Allocate CTA observation
        GCTAObservation obs;

        // Get event file name
        std::string filename = (*this)["infile"].filename();

        // Try first to open as FITS file
        try {

            // Load data
            obs.load(filename);

            // Get other task parameters
            m_stat  = gammalib::toupper((*this)["stat"].string());
            m_caldb = (*this)["caldb"].string();
            m_irf   = (*this)["irf"].string();

            // Set statistics
            obs.statistics(m_stat);

            // Set calibration database. If specified parameter is a
            // directory then use this as the pathname to the calibration
            // database. Otherwise interpret this as the instrument name,
            // the mission being "cta"
            GCaldb caldb;
            if (gammalib::dir_exists(m_caldb)) {
                caldb.rootdir(m_caldb);
            }
            else {
                caldb.open("cta", m_caldb);
            }

            // Set reponse
            obs.response(m_irf, caldb);

            // Append observation to container
            m_obs.append(obs);

        }

        // ... otherwise try to open as XML file
        catch (GException::fits_open_error &e) {

            // Load observations from XML file
            m_obs.load(filename);

            // Check if all observations have response information. If
            // not, get the calibration database parameters and set
            // the response properly
            bool asked_for_response = false;
            for (int i = 0; i < m_obs.size(); ++i) {

                // Get CTA observation
                GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

                // Continue only if observation is a CTA observation
                if (obs != NULL) {

                    // If response is not valid then set response from
                    // task parameters
                    if (!obs->hasresponse()) {

                        // Get calibration parameters
                        if (!asked_for_response) {
                            m_caldb = (*this)["caldb"].string();
                            m_irf   = (*this)["irf"].string();
                            asked_for_response = true;
                        }

                        // Set response
                        GCaldb caldb;
                        if (gammalib::dir_exists(m_caldb)) {
                            caldb.rootdir(m_caldb);
                        }
                        else {
                            caldb.open("cta", m_caldb);
                        }
                        obs->response(m_irf, caldb);

                    } // endif: there was no response
                } // endif: observation was a CTA observation

            } // endfor: looped over observations

        } // endcatch: file was an XML file

    } // endif: there was no observation in the container

    // If there is are no models associated with the observations then
    // load now the model definition
    if (m_obs.models().size() == 0) {

        // Get models XML filename
        std::string filename = (*this)["srcmdl"].filename();

        // Setup models for optimizing.
        m_obs.models(GModels(filename));

    } // endif: no models were associated with observations

    // Get other parameters
    m_refit       = (*this)["refit"].boolean();
    m_tscalc      = (*this)["tscalc"].boolean();
    m_apply_edisp = (*this)["edisp"].boolean();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (m_read_ahead) {
        m_outmdl = (*this)["outmdl"].filename();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Optimize model parameters using Levenberg-Marquardt method
 ***************************************************************************/
void ctlike::optimize_lm(void)
{
    // Free any existing optimizer
    if (m_opt != NULL) delete m_opt;
    m_opt = NULL;

    // Allocate optimizer. The logger is only passed to the optimizer
    // constructor if optimizer logging is requested.
    GOptimizerLM* opt = (logTerse()) ? new GOptimizerLM(log)
                                     : new GOptimizerLM();

    // Assign optimizer
    m_opt = opt;

    // Set optimizer parameters
    opt->max_iter(m_max_iter);
    opt->max_stalls(m_max_stall);

    // Write Header for optimization and indent for optimizer logging
    if (logTerse()) {
        log << std::endl;
        log.header1("Maximum likelihood optimisation");
        log.indent(1);
    }

    // Perform LM optimization
    m_obs.optimize(*opt);

    // Optionally refit
    if (m_refit) {
        m_obs.optimize(*opt);
    }

    // Store maximum log likelihood value
    m_logL = -(opt->value());

    // Write optimization results
    log.indent(0);
    if (logTerse()) {
        log << std::endl;
        log.header1("Maximum likelihood optimization results");
        log << *opt << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Re-optimize model parameters using Levenberg-Marquardt method
 *        for TS computation
 ***************************************************************************/
double ctlike::reoptimize_lm(void)
{
    // Allocate optimizer. The logger is only passed to the optimizer
    // constructor if optimizer logging is requested.
    GOptimizerLM* opt = (logTerse()) ? new GOptimizerLM(log)
                                     : new GOptimizerLM();

    // Set optimizer parameters
    opt->max_iter(m_max_iter);
    opt->max_stalls(m_max_stall);

    // Write Header for optimization and indent for optimizer logging
    if (logTerse()) {
        log << std::endl;
        log.header1("Maximum likelihood re-optimisation");
        log.indent(1);
    }

    // Perform LM optimization
    m_obs.optimize(*opt);

    // Optionally refit
    if (m_refit) {
        m_obs.optimize(*opt);
    }

    // Store maximum log likelihood value
    double logL = -(opt->value());

    // Write optimization results
    log.indent(0);
    if (logTerse()) {
        log << std::endl;
        log.header1("Maximum likelihood re-optimization results");
        log << *opt << std::endl;
    }

    // Return
    return (logL);
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void ctlike::init_members(void)
{
    // Initialise members
    m_stat.clear();
    m_caldb.clear();
    m_irf.clear();
    m_outmdl.clear();
    m_obs.clear();
    m_refit       = false;
    m_tscalc      = false;
    m_max_iter    = 100;   // Set maximum number of iterations
    m_max_stall   = 10;    // Set maximum number of stalls
    m_logL        = 0.0;
    m_opt         = NULL;
    m_read_ahead  = false;
    m_apply_edisp = false;

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
void ctlike::copy_members(const ctlike& app)
{
    // Copy attributes
    m_stat        = app.m_stat;
    m_refit       = app.m_refit;
    m_tscalc      = app.m_tscalc;
    m_caldb       = app.m_caldb;
    m_irf         = app.m_irf;
    m_outmdl      = app.m_outmdl;
    m_obs         = app.m_obs;
    m_max_iter    = app.m_max_iter;
    m_max_stall   = app.m_max_stall;
    m_logL        = app.m_logL;
    m_opt         = app.m_opt->clone();
    m_read_ahead  = app.m_read_ahead;
    m_apply_edisp = app.m_apply_edisp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctlike::free_members(void)
{
    // Free members
    if (m_opt != NULL) delete m_opt;

    // Mark pointers as free
    m_opt = NULL;

    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}
