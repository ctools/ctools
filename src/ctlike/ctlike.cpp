/***************************************************************************
 *                ctlike - Maximum likelihood fitting tool                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @brief Maximum likelihood fitting tool implementation
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
ctlike::ctlike(void) : ctlikelihood(CTLIKE_NAME, CTLIKE_VERSION)
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
 * Constructs ctlike tool from an observations container.
 ***************************************************************************/
ctlike::ctlike(const GObservations& obs) :
        ctlikelihood(CTLIKE_NAME, CTLIKE_VERSION, obs)
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
 *
 * Constructs an instance of the ctlike tool that will parse user parameters
 * that are provided as command line arguments.
 ***************************************************************************/
ctlike::ctlike(int argc, char *argv[]) : 
        ctlikelihood(CTLIKE_NAME, CTLIKE_VERSION, argc, argv)
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
ctlike::ctlike(const ctlike& app) : ctlikelihood(app)
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
 * @return Application.
 ***************************************************************************/
ctlike& ctlike::operator=(const ctlike& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->ctlikelihood::operator=(app);

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
 * @brief Clear ctlike tool
 *
 * Clears ctlike tool.
 ***************************************************************************/
void ctlike::clear(void)
{
    // Free members
    free_members();
    this->ctlikelihood::free_members();
    this->ctobservation::free_members();
    this->ctool::free_members();

    // Clear base class (needed to conserve tool name and version)
    this->GApplication::clear();

    // Initialise members
    this->ctool::init_members();
    this->ctobservation::init_members();
    this->ctlikelihood::init_members();
    init_members();

    // Write header into logger
    log_header();

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

    // Set energy dispersion flags of all CTA observations and save old
    // values in save_edisp vector
    std::vector<bool> save_edisp = set_edisp(m_obs, m_apply_edisp);

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Optimize model parameters using LM optimizer
    optimize_lm();

    // Store copy of curvature matrix
    GMatrixSparse curvature =
        *(const_cast<GObservations::likelihood&>(m_obs.function()).curvature());
    
    // Store Npred
    double npred = m_obs.npred();
    
    // Store models for which TS should be computed
    std::vector<std::string> ts_srcs;
    GModels models_orig = m_obs.models();
    for (int i = 0; i < models_orig.size(); ++i) {
        GModel* model = models_orig[i];
        if (model->tscalc()) {
            ts_srcs.push_back(model->name());
        }
    }

    // Compute TS values if requested
    if (!ts_srcs.empty()) {

        // Store original maximum likelihood and models
        double  logL_src = m_logL;
        GModels models   = m_obs.models();

        // Fix spatial parameters if requested
        if (m_fix_spat_for_ts) {

            // Loop over all models
            for (int i = 0; i < models.size(); ++i) {

                // Continue only if model is skymodel
                GModelSky* sky= dynamic_cast<GModelSky*>(models[i]);
                if (sky != NULL) {

                    // Fix spatial parameters
                    GModelSpatial* spatial = sky->spatial();
                    for (int j = 0; j < spatial->size(); j++) {
                        (*spatial)[j].fix();
                    } // endfor: looped over spatial parameters

                } // endif: there was a sky model

            } // endfor: looped over models

        } // endif: spatial parameter should be fixed

        // Loop over stored models, remove source and refit
        for (int i = 0; i < ts_srcs.size(); ++i) {
            models.remove(ts_srcs[i]);
            m_obs.models(models);    
            double logL_nosrc = reoptimize_lm();
            double ts         = 2.0 * (logL_src-logL_nosrc);
            models_orig[ts_srcs[i]]->ts(ts);
            models = models_orig;
        }

        // Restore best fit values
        m_obs.models(models_orig);
    }

    // Compute number of observed events in all observations
    double num_events = 0.0;
    for (int i = 0; i < m_obs.size(); ++i) {
        double data = m_obs[i]->events()->number();
        if (data >= 0.0) {
            num_events += data;
        }
    }

    // Write results into logger
    log_header1(NORMAL, "Maximum likelihood optimisation results");
    log_string(NORMAL, m_opt.print(m_chatter));
    log_value(NORMAL, "Maximum log likelihood", gammalib::str(m_logL,3));
    log_value(NORMAL, "Observed events  (Nobs)", gammalib::str(num_events,3));
    log_value(NORMAL, "Predicted events (Npred)", gammalib::str(npred,3)+
              " (Nobs - Npred = "+gammalib::str(num_events-npred)+")");
    log_string(NORMAL, m_obs.models().print(m_chatter));

    // Restore energy dispersion flags of all CTA observations
    restore_edisp(m_obs, save_edisp);

    // Restore curvature matrix
    *(const_cast<GObservations::likelihood&>(m_obs.function()).curvature()) =
        curvature;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save results
 *
 * This method saves the fit results into a XML file and a FITS file. If
 * the filename parameters are "NONE", the files are not saved.
 ***************************************************************************/
void ctlike::save(void)
{
    // Write header
    log_header1(TERSE, "Save results");

    // Get output filenames
    m_outmodel  = (*this)["outmodel"].filename();
    m_outcovmat = (*this)["outcovmat"].filename();

    // Save only if filename is valid
    if (is_valid_filename(m_outmodel)) {

        // Log filename
        log_value(NORMAL, "Model definition file", m_outmodel.url());

        // Write results out as XML model
        m_obs.models().save(m_outmodel);

    }

    // ... otherwise signal that file was not saved
    else {
        log_value(NORMAL, "Model definition file", "NONE");
    }

    // Save covariance matrix if filename is valid
    if (is_valid_filename(m_outcovmat)) {

        // Log filename
        log_value(NORMAL, "Covariance matrix file", m_outcovmat.url());

        // Save covarianve matrix
        m_obs.function().save(m_outcovmat);

    }

    // ... otherwise signal that no covariance matrix was not saved
    else {
        log_value(NORMAL, "Covariance matrix file", "NONE");
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
    // Setup observations from "inobs" parameter
    setup_observations(m_obs);

    // If only single observation is used, read statistics parameter
    if (!m_use_xml) {

        // Get other task parameters
        std::string statistics = gammalib::toupper((*this)["stat"].string());

        // Set statistics
        (*m_obs[0]).statistics(statistics);
    }

    // If there is are no models associated with the observations then
    // load now the model definition
    if (m_obs.models().size() == 0) {

        // Get models XML filename
        std::string filename = (*this)["inmodel"].filename();

        // Setup models for optimizing.
        m_obs.models(GModels(filename));

    } // endif: no models were associated with observations

    // Get other parameters
    m_refit           = (*this)["refit"].boolean();
    m_apply_edisp     = (*this)["edisp"].boolean();
    m_fix_spat_for_ts = (*this)["fix_spat_for_ts"].boolean();
    m_chatter         = static_cast<GChatter>((*this)["chatter"].integer());

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outmodel  = (*this)["outmodel"].filename();
        m_outcovmat = (*this)["outcovmat"].filename();
    }

    // Set optimizer logger
    if (logNormal()) {
        static_cast<GOptimizerLM*>(&m_opt)->logger(&log);
    }
    else {
        static_cast<GOptimizerLM*>(&m_opt)->logger(NULL);
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Optimise model parameters
 *
 * Optimise model parameters using a maximum likelihood fit.
 ***************************************************************************/
void ctlike::optimize_lm(void)
{
    // Write header
    log_header1(TERSE, "Maximum likelihood optimisation");
    log.indent(1);

    // Compute number of fitted parameters
    int nfit = 0;
    for (int i = 0; i < m_obs.models().size(); ++i) {
        const GModel* model = m_obs.models()[i];
        for (int k = 0; k < model->size(); ++k) {
            if ((*model)[k].is_free()) {
                nfit++;
            }
        }
    }

    // Notify if all parameters are fixed
    if (nfit == 0) {
        log_string(TERSE, "WARNING: All model parameters are fixed!");
        log_string(TERSE, "         ctlike will proceed without fitting parameters.");
        log_string(TERSE, "         All curvature matrix elements will be zero.");
    }

    // Perform LM optimization
    m_obs.optimize(m_opt);

    // Optionally refit
    if (m_refit) {

        // Dump new header
        log.indent(0);
        log_header1(TERSE, "Maximum likelihood re-optimisation");
        log.indent(1);

        // Optimise again
        m_obs.optimize(m_opt);

    }

    // Optionally show curvature matrix
    log_header1(EXPLICIT, "Curvature matrix");
    log.indent(1);
    log_string(EXPLICIT, (const_cast<GObservations::likelihood&>
                                    (m_obs.function()).curvature())->print());

    // Compute errors
    m_obs.errors(m_opt);

    // Store maximum log likelihood value
    m_logL = -(m_opt.value());

    // Remove indent
    log.indent(0);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Re-optimise model parameters for TS computation
 *
 * Re-optimise the model parameters using a maximum likelihood fit for
 * computation of the Test Statistic value for a given source.
 ***************************************************************************/
double ctlike::reoptimize_lm(void)
{
    // Write Header for optimization and indent for optimizer logging
    log_header1(TERSE, "Maximum likelihood re-optimisation");
    log.indent(1);

    // Create a clone of the optimizer for the re-optimisation
    GOptimizer* opt = m_opt.clone();

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
    log_header1(NORMAL, "Maximum likelihood re-optimisation results");
    log_string(NORMAL, opt->print(m_chatter));

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
    m_outmodel.clear();
    m_outcovmat.clear();
    m_refit           = false;
    m_max_iter        = 100;   // Set maximum number of iterations
    m_max_stall       = 10;    // Set maximum number of stalls
    m_logL            = 0.0;
    m_apply_edisp     = false;
    m_fix_spat_for_ts = false;
    m_chatter         = static_cast<GChatter>(2);

    // Set logger properties
    log.date(true);

    // Set optimizer parameters
    m_opt.max_iter(m_max_iter);
    m_opt.max_stalls(m_max_stall);

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
    m_refit           = app.m_refit;
    m_outmodel        = app.m_outmodel;
    m_outcovmat       = app.m_outcovmat;
    m_max_iter        = app.m_max_iter;
    m_max_stall       = app.m_max_stall;
    m_logL            = app.m_logL;
    m_apply_edisp     = app.m_apply_edisp;
    m_fix_spat_for_ts = app.m_fix_spat_for_ts;
    m_chatter         = app.m_chatter;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctlike::free_members(void)
{
    // Return
    return;
}
