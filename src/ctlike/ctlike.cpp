/***************************************************************************
 *                ctlike - Maximum likelihood fitting tool                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#endif

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
ctlike::ctlike(void) : ctlikelihood(CTLIKE_NAME, VERSION)
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
        ctlikelihood(CTLIKE_NAME, VERSION, obs)
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
        ctlikelihood(CTLIKE_NAME, VERSION, argc, argv)
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

    // Compute number of observed events in all observations
    m_nobs = 0.0;
    for (int i = 0; i < m_obs.size(); ++i) {
        double data = m_obs[i]->nobserved();
        if (data >= 0.0) {
            m_nobs += data;
        }
    }

    // Optimize model parameters using LM optimizer
    optimize_lm();

    // Store copy of curvature matrix
    GMatrixSparse curvature =
        *(const_cast<GObservations::likelihood&>(m_obs.function()).curvature());

    // Store Npred
    m_npred = m_obs.npred();

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

        // Write general fit results in logger (will be repeated at the
        // end)
        log_header1(EXPLICIT, "Maximum likelihood optimisation results");
        log_string(EXPLICIT, m_opt.print(m_chatter));
        log_value(EXPLICIT, "Maximum log likelihood", gammalib::str(m_logL,3));
        log_value(EXPLICIT, "Observed events  (Nobs)", gammalib::str(m_nobs,3));
        log_value(EXPLICIT, "Predicted events (Npred)", gammalib::str(m_npred,3)+
                  " (Nobs - Npred = "+gammalib::str(m_nobs-m_npred)+")");
        log_string(VERBOSE, m_obs.models().print(m_chatter));

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

            // Create copy of models
            GModels models_copy = models;

            // Remove source of interest
            models_copy.remove(ts_srcs[i]);

            // Set models for fitting
            m_obs.models(models_copy);

            // Re-optimise log-likelihood for the source removed
            double logL_nosrc = reoptimize_lm();

            // Compute source TS value
            double ts = 2.0 * (logL_src-logL_nosrc);

            // Store TS value in original model
            models_orig[ts_srcs[i]]->ts(ts);

        } // endfor: looped over sources

        // Restore best fit values
        m_obs.models(models_orig);

    } // endif: requested TS computation

    // Write results into logger
    log_header1(NORMAL, "Maximum likelihood optimisation results");
    log_string(NORMAL, m_opt.print(m_chatter));
    log_value(NORMAL, "Maximum log likelihood", gammalib::str(m_logL,3));
    log_value(NORMAL, "Observed events  (Nobs)", gammalib::str(m_nobs,3));
    log_value(NORMAL, "Predicted events (Npred)", gammalib::str(m_npred,3)+
              " (Nobs - Npred = "+gammalib::str(m_nobs-m_npred)+")");
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
 * This method saves the fit results and the covariance matrix. The fit
 * results are written into a XML file while the covariance matrix is written
 * into either a FITS or a CSV file, depending on the file type extension
 * (an extension of `.fits` or `.fit` produce a FITS file, any other
 * extension produces a CSV file).
 *
 * If the filenames are `NONE` no information is saved.
 ***************************************************************************/
void ctlike::save(void)
{
    // Write header
    log_header1(TERSE, "Save results");

    // Save model only if filename is valid
    if ((*this)["outmodel"].is_valid()) {

        // Generate XML instance
        GXml xml = xml_result();

        // Get output filename
        m_outmodel = (*this)["outmodel"].filename();

        // Log filename
        log_value(NORMAL, "Model definition file", m_outmodel.url());

        // Write results out as XML model
        xml.save(m_outmodel);

    } // endif: filename was valid

    // ... otherwise signal that file was not saved
    else {
        log_value(NORMAL, "Model definition file", "NONE");
    }

    // Save covariance matrix only if filename is valid
    if ((*this)["outcovmat"].is_valid()) {

        // Get output filenames
        m_outcovmat = (*this)["outcovmat"].filename();

        // Log filename
        log_value(NORMAL, "Covariance matrix file", m_outcovmat.url());

        // Save covariance matrix
        m_obs.function().save(m_outcovmat);

    }

    // ... otherwise signal that no covariance matrix was not saved
    else {
        log_value(NORMAL, "Covariance matrix file", "NONE");
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
void ctlike::init_members(void)
{
    // Initialise members
    m_outmodel.clear();
    m_outcovmat.clear();
    m_refit           = false;
    m_apply_edisp     = false;
    m_fix_spat_for_ts = false;
    m_chatter         = static_cast<GChatter>(2);
    m_logL            = 0.0;
    m_nobs            = 0.0;
    m_npred           = 0.0;

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
    m_refit           = app.m_refit;
    m_outmodel        = app.m_outmodel;
    m_outcovmat       = app.m_outcovmat;
    m_apply_edisp     = app.m_apply_edisp;
    m_fix_spat_for_ts = app.m_fix_spat_for_ts;
    m_chatter         = app.m_chatter;
    m_logL            = app.m_logL;
    m_nobs            = app.m_nobs;
    m_npred           = app.m_npred;

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


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all required task parameters from the parameter file.
 ***************************************************************************/
void ctlike::get_parameters(void)
{
    // Setup observations from "inobs" parameter
    setup_observations(m_obs);

    // Set observation statistic
    set_obs_statistic(gammalib::toupper((*this)["statistic"].string()));

    // If there is are no models associated with the observations then
    // load now the model definition
    if (m_obs.models().size() == 0) {

        // Get models XML filename
        std::string filename = (*this)["inmodel"].filename();

        // Setup models for optimizing.
        m_obs.models(GModels(filename));

    } // endif: no models were associated with observations

    // Set optimizer characteristics from user parameters
    m_opt.eps((*this)["like_accuracy"].real());
    m_opt.max_iter((*this)["max_iter"].integer());

    // Get other parameters
    m_refit           = (*this)["refit"].boolean();
    m_apply_edisp     = (*this)["edisp"].boolean();
    m_fix_spat_for_ts = (*this)["fix_spat_for_ts"].boolean();
    m_chatter         = static_cast<GChatter>((*this)["chatter"].integer());

    // If needed later, query output filenames now
    if (read_ahead()) {
        (*this)["outmodel"].query();
        (*this)["outcovmat"].query();
    }

    // Set optimizer logger
    if (logNormal()) {
        static_cast<GOptimizerLM*>(&m_opt)->logger(&log);
    }
    else {
        static_cast<GOptimizerLM*>(&m_opt)->logger(NULL);
    }

    // Set number of OpenMP threads
    #ifdef _OPENMP
    int nthreads = (*this)["nthreads"].integer();
    if (nthreads > 0) {
        omp_set_num_threads(nthreads);
    }
    #endif

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
    log_header1(EXPLICIT, "Maximum likelihood re-optimisation results");
    log_string(EXPLICIT, opt->print(m_chatter));
    log_value(EXPLICIT, "Maximum log likelihood", gammalib::str(logL,3));
    log_value(EXPLICIT, "Observed events  (Nobs)", gammalib::str(m_nobs,3));
    log_value(EXPLICIT, "Predicted events (Npred)", gammalib::str(m_obs.npred(),3)+
              " (Nobs - Npred = "+gammalib::str(m_nobs-m_obs.npred())+")");
    log_string(VERBOSE, m_obs.models().print(m_chatter));

    // Return
    return (logL);
}


/***********************************************************************//**
 * @brief Generate XML result
 *
 * @return XML result
 *
 * Generates the XML result composed of the ctlike results and the model
 * fitting results.
 ***************************************************************************/
GXml ctlike::xml_result(void) const
{
    // Initialise XML result
    GXml xml;

    // Set fit status
    std::string status;
    switch (m_opt.status()) {
    case G_LM_CONVERGED:
        status.append("converged");
        break;
    case G_LM_STALLED:
        status.append("stalled");
        break;
    case G_LM_SINGULAR:
        status.append("singular curvature matrix encountered");
        break;
    case G_LM_NOT_POSTIVE_DEFINITE:
        status.append("curvature matrix not positive definite");
        break;
    case G_LM_BAD_ERRORS:
        status.append("errors are inaccurate");
        break;
    default:
        status.append("unknown");
        break;
    }

    // Set flag strings
    std::string refit           = (m_refit) ? "yes" : "no";
    std::string edisp     = (m_apply_edisp) ? "yes" : "no";
    std::string fix_spat_for_ts = (m_fix_spat_for_ts) ? "yes" : "no";

    // Write ctlike results into XML instance
    if (xml.elements("ctlike_results") == 0) {
        xml.append(GXmlElement("ctlike_results title=\"ctlike fit results\""));
    }
    GXmlElement* result = xml.element("ctlike_results", 0);
    result->append(GXmlElement("status", status));
    result->append(GXmlElement("log-likelihood", m_logL));
    result->append(GXmlElement("precision", m_opt.eps()));
    result->append(GXmlElement("iterations", m_opt.iter()));
    result->append(GXmlElement("lambda", m_opt.lambda()));
    result->append(GXmlElement("total_parameters", m_opt.npars()));
    result->append(GXmlElement("fitted_parameters", m_opt.nfree()));
    result->append(GXmlElement("observed-events", m_nobs));
    result->append(GXmlElement("predicted-events", m_npred));
    result->append(GXmlElement("refit", refit));
    result->append(GXmlElement("edisp", edisp));
    result->append(GXmlElement("fix-spatial-for-ts", fix_spat_for_ts));
    result->append(GXmlElement("elapsed-time", this->telapse()));
    result->append(GXmlElement("cpu-seconds", this->celapse()));

    // Write model results into XML instance
    m_obs.models().write(xml);

    // Return XML result
    return xml;
}
