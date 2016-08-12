/***************************************************************************
 *                 cterror - Parameter error calculation tool              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2016 by Florent Forest                              *
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
 * @file cterror.cpp
 * @brief Parameter error calculation tool interface implementation
 * @author Florent Forest
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "cterror.hpp"
#include "GTools.hpp"
#include "GOptimizer.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ERR_BISECTION         "cterror::error_bisection(double&, double&)"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty cterror tool.
 ***************************************************************************/
cterror::cterror(void) : ctlikelihood(CTERROR_NAME, CTERROR_VERSION)
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
 * Constructs cterror tool from an observations container.
 ***************************************************************************/
cterror::cterror(const GObservations& obs) :
         ctlikelihood(CTERROR_NAME, CTERROR_VERSION, obs)
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
 * Constructs an instance of the cterror tool that will parse user
 * parameters that are provided as command line arguments.
 ***************************************************************************/
cterror::cterror(int argc, char *argv[]) :
         ctlikelihood(CTERROR_NAME, CTERROR_VERSION, argc, argv)
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
 *
 * Constructs an instance of the cterror tool by copying information from
 * another ctulimit tool.
 ***************************************************************************/
cterror::cterror(const cterror& app) : ctlikelihood(app)
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
 *
 * Destructs the cterror tool.
 ***************************************************************************/
cterror::~cterror(void)
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
 *
 * Assigns a cterror tool.
 ***************************************************************************/
cterror& cterror::operator=(const cterror& app)
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
 * @brief Clear cterror tool
 *
 * Clears cterror tool.
 ***************************************************************************/
void cterror::clear(void)
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
 * @brief Compute parameter errors using a likelihood profile method
 *
 * Computes the parameter errors using a likelihood profile method.
 ***************************************************************************/
void cterror::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Set energy dispersion flags of all CTA observations and save old
    // values in save_edisp vector
    std::vector<bool> save_edisp = set_edisp(m_obs, m_apply_edisp);

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write header into logger
    log_header1(TERSE, "Compute best-fit likelihood");

    // Optimize and save best log-likelihood
    m_obs.optimize(m_opt);
    m_obs.errors(m_opt);
    m_best_logL = m_obs.logL();

    // Store optimizer for later recovery
    GOptimizerLM best_opt = m_opt;

    // Write optimisation results and models into logger
    log_string(NORMAL, m_opt.print(m_chatter));
    log_value(NORMAL, "Maximum log likelihood", gammalib::str(m_best_logL,3));
    log_string(NORMAL, m_obs.models().print(m_chatter));

    // Save best fitting models
    GModels models_best = m_obs.models();

    // Get pointer on model
    GModel* model = models_best[m_srcname];

    // Get number of parameters
    int npars = model->size();

    // Loop over parameters of sky model
    for (int i = 0; i < npars; ++i) {

        // Skip parameter if it is fixed
        if (model->at(i).is_fixed()) {
            continue;
        }

        // Initialise with best fitting models
        m_obs.models(models_best);

        // Get pointer on model parameter
        GModels& current_models = const_cast<GModels&>(m_obs.models());
        m_model_par             = &(current_models[m_srcname]->at(i));

        // Extract current value
        m_value = m_model_par->factor_value();

        // Compute parameter bracketing
        double parmin = std::max(m_model_par->factor_min(),
                                 m_value - 10.0*m_model_par->factor_error());
        double parmax = std::min(m_model_par->factor_max(),
                                 m_value + 10.0*m_model_par->factor_error());

        // Write header and initial parameters into logger
        log_header1(TERSE, "Compute error for source \""+m_srcname+"\""
                           " parameter \""+m_model_par->name()+"\"");
        log_value(NORMAL, "Confidence level",
                  gammalib::str(m_confidence*100.0)+" %");
        log_value(NORMAL, "Log-likelihood difference", m_dlogL);
        log_value(NORMAL, "Initial factor range",
                  "["+gammalib::str(parmin)+", "+gammalib::str(parmax)+"]");

        // Compute lower and upper boundaries
        double value_lo = error_bisection(parmin, m_value);
        double value_hi = error_bisection(m_value, parmax);

        // Compute errors
        double error           = 0.5 * (value_hi - value_lo);
        double error_neg       = m_value  - value_lo;
        double error_pos       = value_hi - m_value;
        double error_value     = std::abs(error*m_model_par->scale());
        double error_value_neg = std::abs(error_neg*m_model_par->scale());
        double error_value_pos = std::abs(error_pos*m_model_par->scale());

        // Write results into logger
        std::string unit = " " + m_model_par->unit();
        log_value(NORMAL, "Lower parameter factor", value_lo);
        log_value(NORMAL, "Upper parameter factor", value_hi);
        log_value(NORMAL, "Error from curvature",
                  gammalib::str(m_model_par->error()) + unit);
        log_value(NORMAL, "Error from profile",
                  gammalib::str(error_value) + unit);
        log_value(NORMAL, "Negative profile error",
                  gammalib::str(error_value_neg) + unit);
        log_value(NORMAL, "Positive profile error",
                  gammalib::str(error_value_pos) + unit);

        // Save error result
        model->at(i).factor_error(error);

    } // endfor: looped over spectral parameters

    // Restore best fitting models (now with new errors computed)
    m_obs.models(models_best);

    // Recover optimizer
    m_opt = best_opt;

    // Restore energy dispersion flags of all CTA observations
    restore_edisp(m_obs, save_edisp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save model
 *
 * Saves the model into an XML file.
 ***************************************************************************/
void cterror::save(void)
{
    // Write header
    log_header1(TERSE, "Save results");

    // Get output filename
    m_outmodel = (*this)["outmodel"].filename();

    // Save only if filename is valid
    if (is_valid_filename(m_outmodel)) {

        // Log filename
        log_value(NORMAL, "Model definition file", m_outmodel.url());

        // Write results out as XML model
        m_obs.models().save(m_outmodel.url());

    }

    // ... otherwise signal that file was not saved
    else {
        log_value(NORMAL, "Model definition file", "NONE");
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
void cterror::init_members(void)
{
    // Initialise user parameters
    m_srcname.clear();
    m_outmodel.clear();
    m_confidence  = 0.68;
    m_tol         = 1.0e-3;
    m_max_iter    = 50;
    m_apply_edisp = false;
    m_chatter     = static_cast<GChatter>(2);

    // Initialise protected members
    m_value       = 0.0;
    m_dlogL       = 0.0;
    m_best_logL   = 0.0;
    m_model_par   = NULL;

    // Set optimizer parameters
    m_opt.max_iter(m_max_iter);
    m_opt.max_stalls(10);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void cterror::copy_members(const cterror& app)
{
    // Copy user parameters
    m_srcname     = app.m_srcname;
    m_outmodel    = app.m_outmodel;
    m_confidence  = app.m_confidence;
    m_tol         = app.m_tol;
    m_max_iter    = app.m_max_iter;
    m_apply_edisp = app.m_apply_edisp;
    m_chatter     = app.m_chatter;

    // Copy protected members
    m_value     = app.m_value;
    m_dlogL     = app.m_dlogL;
    m_best_logL = app.m_best_logL;
    m_model_par = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void cterror::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user.
 ***************************************************************************/
void cterror::get_parameters(void)
{
    // Setup observations from "inobs" parameter
    setup_observations(m_obs);
    
    // Setup models from "inmodel" parameter
    setup_models(m_obs, (*this)["srcname"].string());

    // Get name of test source
    m_srcname = (*this)["srcname"].string();

    // Read energy dispersion flag
    m_apply_edisp = (*this)["edisp"].boolean();

    // Get confidence level and transform into log-likelihood difference
    m_confidence = (*this)["confidence"].real();
    double sigma = gammalib::erfinv(m_confidence) * gammalib::sqrt_two;
    m_dlogL      = (sigma*sigma) / 2.0;

    // Read other parameters
    m_tol      = (*this)["tol"].real();
    m_max_iter = (*this)["max_iter"].integer();
    m_chatter  = static_cast<GChatter>((*this)["chatter"].integer());

    // Read ahead parameters that are only needed when the tool gets
    // executed
    if (read_ahead()) {
        m_outmodel = (*this)["outmodel"].filename();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Calculate error using a bisection method
 *
 * @param[in] min Minimum parameter value
 * @param[in] max Maximum parameter value
 *
 * Calculates the error using a bisection method.
 ***************************************************************************/
double cterror::error_bisection(const double& min, const double& max)
{
    // Copy values to working values
    double wrk_min = min;
    double wrk_max = max;

    // Initialise iteration counter
    int iter = 1;

    // Initialize mid value
    double mid = (wrk_min + wrk_max) / 2.0;

    // Loop until breaking condition is reached
    while (true) {

        // Throw exception if maximum iterations are reached
        if (iter > m_max_iter) {
            if (wrk_min - m_model_par->factor_min() < m_tol) {
                std::string msg = "The \""+m_model_par->name()+"\" parameter "
                                  "minimum has been reached during error "
                                  "calculation. To obtain accurate errors, "
                                  "consider setting the minimum parameter "
                                  "value to a lower value, and re-run "
                                  "cterror.";
                log_string(TERSE, msg);
                break;
            }
            else if (m_model_par->factor_max() - wrk_max < m_tol) {
                std::string msg = "The \""+m_model_par->name()+"\" parameter "
                                  "maximum has been reached during error "
                                  "calculation. To obtain accurate errors, "
                                  "consider setting the maximum parameter "
                                  "value to a higher value, and re-run "
                                  "cterror.";
                log_string(TERSE, msg);
                break;
            }
            else {
                std::string msg = "The maximum number of "+
                                  gammalib::str(m_max_iter)+" iterations has "
                                  "been reached. Please increase the "
                                  "\"max_iter\" parameter, and re-run "
                                  "cterror.";
                throw GException::invalid_value(G_ERR_BISECTION, msg);
            }
        }

        // Compute center of boundary
        mid = (wrk_min + wrk_max) / 2.0;

        // Calculate function value
        double eval_mid = evaluate(*m_model_par, mid) - (m_best_logL + m_dlogL);

        // Write interval into logger
        log_value(EXPLICIT, "  Iteration "+gammalib::str(iter),
                  "["+gammalib::str(wrk_min)+", "+gammalib::str(wrk_max)+"]");

        // Check for convergence inside tolerance
        if (std::abs(eval_mid) < m_tol) {
            break;
        }

        // Check if interval is smaller than 1.0e-6
        if (std::abs(wrk_max-wrk_min) < 1.0e-6) {
            break;
        }

        // If we are on the crescent side of the parabola ...
        if (mid > m_value) {

            // Change boundaries for further iteration
            if (eval_mid > 0.0) {
                wrk_max = mid;
            }
            else if (eval_mid < 0.0) {
                wrk_min = mid;
            }
        }

        // ... otherwise we are on the decrescent side of the parabola
        else {

            // Change boundaries for further iteration
            if (eval_mid > 0.0) {
                wrk_min = mid;
            }
            else if (eval_mid < 0.0) {
                wrk_max = mid;
            }
        }

        // Increment counter
        iter++;

    } // endwhile

    // Return mid value
    return mid;

}
