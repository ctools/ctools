/***************************************************************************
 *                 cterror - Parameter error calculation tool              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Florent Forest                                   *
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
#define G_GET_PARAMETERS                         "cterror::get_parameters()"
#define G_ERR_BISECTION         "cterror::error_bisection(double&, double&)"
#define G_EVALUATE                              "cterror::evaluate(double&)"

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
cterror::cterror(void) : ctool(CTERROR_NAME, CTERROR_VERSION)
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
 * This method creates an instance of the class by copying an existing
 * observations container.
 ***************************************************************************/
cterror::cterror(const GObservations& obs) :
         ctool(CTERROR_NAME, CTERROR_VERSION)
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
cterror::cterror(int argc, char *argv[]) :
         ctool(CTERROR_NAME, CTERROR_VERSION, argc, argv)
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
cterror::cterror(const cterror& app) : ctool(app)
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
 ***************************************************************************/
cterror& cterror::operator=(const cterror& app)
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
void cterror::clear(void)
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
 * @brief Computes 
 ***************************************************************************/
void cterror::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
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

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Compute best-fit likelihood");
    }

    // Optimize and save best log-likelihood
    m_obs.optimize(m_opt);
    m_obs.errors(m_opt);
    m_best_logL = m_obs.logL();

    // Store optimizer for later recovery
    GOptimizerLM best_opt = m_opt;

    // Write optimised model into logger
    if (logTerse()) {
        log << m_opt << std::endl;
        log << gammalib::parformat("Maximum log likelihood");
        log << gammalib::str(m_best_logL,3) << std::endl;
        log << m_obs.models() << std::endl;
    }

    // Continue only if source model exists
    if (m_obs.models().contains(m_srcname)) {

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

            // Write header
            if (logTerse()) {
                log << std::endl;
                log.header1("Compute error for source \""+m_srcname+"\""
                            " parameter \""+m_model_par->name()+"\"");
                log << gammalib::parformat("Confidence level");
                log << m_confidence*100.0 << "%" << std::endl;
                log << gammalib::parformat("Log-likelihood difference");
                log << m_dlogL << std::endl;
                log << gammalib::parformat("Initial factor range");
                log << "[";
                log << parmin;
                log << ", ";
                log << parmax;
                log << "]" << std::endl;
            }

            // Compute lower boundary
            double value_lo = error_bisection(parmin, m_value);

            // Write lower parameter value
            if (logTerse()) {
                log << gammalib::parformat("Lower parameter factor");
                log << value_lo << std::endl;
            }

            // Compute upper boundary
            double value_hi = error_bisection(m_value, parmax);

            // Write upper parameter value
            if (logTerse()) {
                log << gammalib::parformat("Upper parameter factor");
                log << value_hi << std::endl;
            }

            // Compute errors
            double error     = 0.5 * (value_hi - value_lo);
            double error_neg = m_value  - value_lo;
            double error_pos = value_hi - m_value;
            //double error_max = std::max(value_hi-m_value, m_value-value_lo);
            //double error_min = std::min(value_hi-m_value, m_value-value_lo);

            // Write errors
            if (logTerse()) {
                log << gammalib::parformat("Error from curvature");
                log << m_model_par->error();
                log << " " << m_model_par->unit() << std::endl;
                log << gammalib::parformat("Error from profile");
                log << std::abs(error*m_model_par->scale());
                log << " " << m_model_par->unit() << std::endl;
                log << gammalib::parformat("Negative profile error");
                log << std::abs(error_neg*m_model_par->scale());
                log << " " << m_model_par->unit() << std::endl;
                log << gammalib::parformat("Positive profile error");
                log << std::abs(error_pos*m_model_par->scale());
                log << " " << m_model_par->unit() << std::endl;
            }

            // Save error result
            model->at(i).factor_error(error);

        } // endfor: looped over spectral parameters

        // Restore best fitting models (now with new errors computed)
        m_obs.models(models_best);

    } // endif: source model exists

    // Recover optimizer
    m_opt = best_opt;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save source model
 *
 * This method saves the soutce model to an XML file.
 ***************************************************************************/
void cterror::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save results");
    }

    // Get output filename
    std::string outmodel = (*this)["outmodel"].filename();

     // Write results out as XML model
    if (gammalib::toupper(outmodel) != "NONE") {
        m_obs.models().save(outmodel);
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
    m_confidence = 0.68;
    m_tol        = 1.0e-3;
    m_max_iter   = 50;
    m_value      = 0.0;

    // Initialise protected members
    m_obs.clear();
    m_opt.clear();
    m_dlogL        = 0.0;
    m_best_logL    = 0.0;
    m_model_par    = NULL;

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
    m_srcname    = app.m_srcname;
    m_confidence = app.m_confidence;
    m_tol        = app.m_tol;
    m_max_iter   = app.m_max_iter;

    // Copy protected members
    m_obs       = app.m_obs;
    m_dlogL     = app.m_dlogL;
    m_best_logL = app.m_best_logL;
    m_opt       = app.m_opt;
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
 * @exception GException::invalid_value
 *            Test source not found.
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user.
 ***************************************************************************/
void cterror::get_parameters(void)
{
    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {

        // Throw exception if no input observation file is given
        require_inobs(G_GET_PARAMETERS);

        // Build observation container
        m_obs = get_observations();

    } // endif: there was no observation in the container


    // If there are no models associated with the observations then
    // load now the model definition
    if (m_obs.models().size() == 0) {

        // Get models XML filename
        std::string filename = (*this)["inmodel"].filename();

        // Setup models for optimizing.
        m_obs.models(GModels(filename));

    } // endif: no models were associated with observations

    // Get name of test source and check container for this name
    m_srcname = (*this)["srcname"].string();
    if (!m_obs.models().contains(m_srcname)) {
        std::string msg = "Source \""+m_srcname+"\" not found in model "
                          "container. Please add a source with that name "
                          "or check for possible typos.";
        throw GException::invalid_value(G_GET_PARAMETERS, msg);
    }

    // Get confidence level and transform into log-likelihood difference
    m_confidence = (*this)["confidence"].real();
    double sigma = gammalib::erfinv(m_confidence) * gammalib::sqrt_two;
    m_dlogL      = (sigma*sigma) / 2.0;

    // Read other parameters
    m_tol      = (*this)["tol"].real();
    m_max_iter = (*this)["max_iter"].integer();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Performs error computation by using a bisection method
 *
 * @param[in] min Minimum parameter value
 * @param[in] max Maximum parameter value
 *
 * This method calculates the error using a bisection method.
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
                if (logTerse()) {
                    log << msg;
                }
                break;
            }
            else if (m_model_par->factor_max() - wrk_max < m_tol) {
                std::string msg = "The \""+m_model_par->name()+"\" parameter "
                                  "maximum has been reached during error "
                                  "calculation. To obtain accurate errors, "
                                  "consider setting the maximum parameter "
                                  "value to a higher value, and re-run "
                                  "cterror.";
                if (logTerse()) {
                    log << msg;
                }
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
        double eval_mid = evaluate(mid);

        // Log interval
        if (logExplicit()) {
            log << gammalib::parformat("  Iteration "+gammalib::str(iter));
            log << "[";
            log << wrk_min;
            log << ", ";
            log << wrk_max;
            log << "] ";
            log << eval_mid;
            log << std::endl;
        }

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


/***********************************************************************//**
 * @brief Evaluates the log-likelihood
 *
 * @param[in] value Parameter factor value
 * @return Log-likelihood value
 *
 * This method evaluates the log-likelihood as a function of the parameter
 * of interest.
 ***************************************************************************/
double cterror::evaluate(const double& value)
{
    // Initialise log-likelihood value
    double logL = 0.0;

    // Check if given parameter is within boundaries
    if (value > m_model_par->factor_min() && value < m_model_par->factor_max()) {

        // Change parameter factor
        m_model_par->factor_value(value);

        // Fix parameter
        m_model_par->fix();

        // Re-optimize
        m_obs.optimize(m_opt);

        // Retrieve likelihood
        logL = m_obs.logL();

        // Free parameter
        m_model_par->free();

    } // endif: value was inside allowed range

    // ... otherwise signal that the parameter went outside the boundaries
    else {
        std::string msg = "Value of parameter \""+m_model_par->name()+"\" "
                          "outside of validity range requested. To omit "
                          "this error please enlarge the parameter range "
                          "in the model XML file.";
        throw GException::invalid_value(G_EVALUATE, msg);
    }

    // Compute function value
    double logL_difference = logL - m_best_logL - m_dlogL;

    // Return
    return logL_difference;
}
