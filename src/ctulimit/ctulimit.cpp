/***************************************************************************
 *                   ctulimit - Upper limit calculation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2022 by Michael Mayer                               *
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
 * @file ctulimit.hpp
 * @brief Upper limit calculation tool interface implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctulimit.hpp"
#include "GTools.hpp"
#include "GOptimizer.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_MODEL_PARAMETER               "ctulimit::get_model_parameter()"
#define G_GET_PARAMETER_BRACKETS "ctulimit::get_parameter_brackets(double&, "\
                                                                   "double&)"
#define G_UL_BISECTION             "ctulimit::ul_bisection(double&, double&)"

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
 * Constructs an empty ctulimit tool.
 ***************************************************************************/
ctulimit::ctulimit(void) : ctlikelihood(CTULIMIT_NAME, VERSION)
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
 * Constructs ctulimit tool from an observations container.
 ***************************************************************************/
ctulimit::ctulimit(const GObservations& obs) :
          ctlikelihood(CTULIMIT_NAME, VERSION, obs)
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
 * Constructs an instance of the ctulimit tool that will parse user
 * parameters that are provided as command line arguments.
 ***************************************************************************/
ctulimit::ctulimit(int argc, char *argv[]) :
          ctlikelihood(CTULIMIT_NAME, VERSION, argc, argv)
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
 * Constructs an instance of the ctulimit tool by copying information from
 * another ctulimit tool.
 ***************************************************************************/
ctulimit::ctulimit(const ctulimit& app) : ctlikelihood(app)
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
 * Destructs the ctulimit tool.
 ***************************************************************************/
ctulimit::~ctulimit(void)
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
 * Assigns a ctulimit tool.
 ***************************************************************************/
ctulimit& ctulimit::operator=(const ctulimit& app)
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
 * @brief Clear ctulimit tool
 *
 * Resets the ctulimit tool to a clean initial state.
 ***************************************************************************/
void ctulimit::clear(void)
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
 * @brief Compute upper limit
 *
 * Computes the upper limit depending on the given algorithm
 ***************************************************************************/
void ctulimit::process(void)
{
    // Get task parameters
    get_parameters();

    // Set energy dispersion flags of all CTA observations and save old
    // values in save_edisp vector
    std::vector<bool> save_edisp = set_edisp(m_obs, m_apply_edisp);

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Save original models
    GModels models_orig = m_obs.models();

    // Initialise best fit optimizer
    GOptimizerLM best_opt;

    // Save original log-likelihood. If the value is zero it has never been
    // computed hence we compute it now. 
    m_best_logL = m_obs.logL();
    if (m_best_logL == 0.0) {

        // Write header into logger
        log_header1(TERSE, "Compute best-fit likelihood");

        // Make sure that requested model parameters is free
        m_model_par->free();

        // Optimize and save best log-likelihood
        m_obs.optimize(m_opt);
        m_obs.errors(m_opt);
        m_best_logL = m_obs.logL();

        // Store optimizer for later recovery
        best_opt = m_opt;

        // Write optimised model into logger
        log_string(NORMAL, m_opt.print(m_chatter));
        log_value(NORMAL,"Maximum log likelihood",gammalib::str(m_best_logL,3));
        log_string(NORMAL, m_obs.models().print(m_chatter));

    } // endif: likelihood was zero

    // Get parameter brackets. Note that this method recomputes m_best_logL
    // to make sure that the value is coherent for the upper limit search.
    double parmin;
    double parmax;
    get_parameter_brackets(parmin, parmax);

    // Write header into logger
    log_header1(TERSE, "Compute upper limit");

    // Write some parameters into logger
    log_value(NORMAL, "Model name", m_skymodel->name());
    log_value(NORMAL, "Parameter name", m_model_par->name());
    log_value(NORMAL, "Parameter value", m_model_par->value());
    log_value(NORMAL, "Confidence level", gammalib::str(m_confidence*100.0)+"%");
    log_value(NORMAL, "Log-likelihood difference", m_dlogL);
    log_value(NORMAL, "Initial factor value range",
              "["+gammalib::str(parmin)+", "+gammalib::str(parmax)+"]");

    // Compute upper limit
    ulimit_bisection(parmin, parmax);

    // Write final parameter into logger
    log_value(NORMAL, "Upper limit factor value", m_model_par->factor_value());
    log_value(NORMAL, "Upper limit param. value", m_model_par->value());

    // Compute upper limit intensity and fluxes
    compute_ulimit();

    // Write header into logger
    log_header1(TERSE, "Upper limit results");

    // Write result into logger
    log_value(TERSE, "Differential flux limit",
              gammalib::str(m_diff_ulimit)+" ph/cm2/s/MeV at "+
              gammalib::str(m_eref)+" TeV");
    log_value(TERSE, "Integral flux limit",
              gammalib::str(m_flux_ulimit)+" ph/cm2/s within ["+
              gammalib::str(m_emin)+"-"+
              gammalib::str(m_emax)+"] TeV");
    log_value(TERSE, "Energy flux limit",
              gammalib::str(m_eflux_ulimit)+" erg/cm2/s within ["+
              gammalib::str(m_emin)+"-"+
              gammalib::str(m_emax)+"] TeV");

    // Recover original models
    m_obs.models(models_orig);

    // Recover optimizer
    m_opt = best_opt;

    // Restore energy dispersion flags of all CTA observations
    restore_edisp(m_obs, save_edisp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save upper limits
 *
 * Saves the upper limit to an ASCII file in comma-separated value format.
 *
 * @todo No yet implemented as we have no clear use case yet for saving
 * the result.
 ***************************************************************************/
void ctulimit::save(void)
{
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
void ctulimit::init_members(void)
{
    // Initialise user parameters
    m_srcname.clear();
    m_parname.clear();
    m_confidence  = 0.95;
    m_sigma_min   = 0.0;
    m_sigma_max   = 0.0;
    m_eref        = 0.0;
    m_emin        = 0.0;
    m_emax        = 0.0;
    m_tol         = 1.0e-6;
    m_max_iter    = 50;
    m_apply_edisp = false;
    m_chatter     = static_cast<GChatter>(2);

    // Initialise protected members
    m_dlogL        = 0.0;
    m_skymodel     = NULL;
    m_model_par    = NULL;
    m_best_logL    = 0.0;
    m_best_value   = 0.0;
    m_best_error   = 0.0;
    m_flux_ulimit  = 0.0;
    m_diff_ulimit  = 0.0;
    m_eflux_ulimit = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctulimit::copy_members(const ctulimit& app)
{
    // Copy user parameters
    m_srcname     = app.m_srcname;
    m_parname     = app.m_parname;
    m_confidence  = app.m_confidence;
    m_sigma_min   = app.m_sigma_min;
    m_sigma_max   = app.m_sigma_max;
    m_eref        = app.m_eref;
    m_emin        = app.m_emin;
    m_emax        = app.m_emax;
    m_tol         = app.m_tol;
    m_max_iter    = app.m_max_iter;
    m_apply_edisp = app.m_apply_edisp;
    m_chatter     = app.m_chatter;

    // Copy protected members
    m_dlogL        = app.m_dlogL;
    m_best_logL    = app.m_best_logL;
    m_best_value   = app.m_best_value;
    m_best_error   = app.m_best_error;
    m_diff_ulimit  = app.m_diff_ulimit;
    m_flux_ulimit  = app.m_flux_ulimit;
    m_eflux_ulimit = app.m_eflux_ulimit;

    // Extract model parameter
    get_model_parameter();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctulimit::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file.
 ***************************************************************************/
void ctulimit::get_parameters(void)
{
    // Setup observations from "inobs" parameter
    setup_observations(m_obs);

    // Set observation statistic
    set_obs_statistic(gammalib::toupper((*this)["statistic"].string()));

    // Setup models from "inmodel" parameter
    setup_models(m_obs, (*this)["srcname"].string());

    // Get name of test source and optional parameter
    m_srcname = (*this)["srcname"].string();
    m_parname = (*this)["parname"].string();

    // Get relevant model and parameter for upper limit computation
    get_model_parameter();

    // Read energy dispersion flag
    m_apply_edisp = (*this)["edisp"].boolean();

    // Get confidence level and transform into log-likelihood difference
    m_confidence = (*this)["confidence"].real();
    double sigma = gammalib::erfinv(m_confidence) * gammalib::sqrt_two;
    m_dlogL      = (sigma*sigma) / 2.0;

    // Read starting boundaries for bisection
    m_sigma_min = (*this)["sigma_min"].real();
    m_sigma_max = (*this)["sigma_max"].real();

    // Read energy values
    m_eref = (*this)["eref"].real();
    m_emin = (*this)["emin"].real();
    m_emax = (*this)["emax"].real();

    // Set optimizer characteristics from user parameters
    m_opt.eps((*this)["like_accuracy"].real());
    m_opt.max_iter((*this)["max_iter"].integer());

    // Read precision
    m_tol      = (*this)["tol"].real();
    m_max_iter = (*this)["max_iter"].integer();

    // Get remaining parameters
    m_chatter = static_cast<GChatter>((*this)["chatter"].integer());

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * @exception GException::invalid_value
 *            Did not find a valid model or model parameter
 *
 * Extracts a pointer to the sky model (m_skymodel) and a pointer to the
 * relevant model parameter (m_model_par) from the model container. If no
 ***************************************************************************/
void ctulimit::get_model_parameter(void)
{
    // Define list of valid spectral model parameters, terminated with a
    // NULL pointer to signal the end of the list
    const char* pars[] = {"Normalization", "Value", "Prefactor",
                          "Integral", "PhotonFlux", "EnergyFlux",
                          NULL};

    // Initialises sky model and model parameters pointers
    m_skymodel  = NULL;
    m_model_par = NULL;

    // If the model container does not container the specified source then
    // throw an exception
    if (!m_obs.models().contains(m_srcname)) {
        std::string msg = "Source \""+m_srcname+"\" not found in model "
                          "container. Please add a source with that name "
                          "or check for possible typos.";
        throw GException::invalid_value(G_GET_MODEL_PARAMETER, msg);
    }

    // Get relevant sky model for upper limit computation. If the model is not
    // a sky model or if the model has a "NodeFunction" as spectral component
    // then throw an exception.
    GModels& models = const_cast<GModels&>(m_obs.models());
    m_skymodel      = dynamic_cast<GModelSky*>(models[m_srcname]);
    if (m_skymodel == NULL) {
        std::string msg = "Source \""+m_srcname+"\" is not a sky model. Please "
                          "specify the name of a sky model for upper limit "
                          "computation.";
        throw GException::invalid_value(G_GET_MODEL_PARAMETER, msg);
    }

    // If a parameter name was specified then use that parameter
    if (!m_parname.empty()) {
        if (m_skymodel->spectral()->has_par(m_parname)) {
            m_model_par = &(m_skymodel->spectral()->operator[](m_parname));
        }
        else {
            std::string msg = "Spectral model of source \""+m_srcname+"\" has "
                              "no parameter \""+m_parname+"\". Please specify "
                              "a value parameter name of leave the parameter "
                              "name blank for autodetermination of an "
                              "intensity- or flux-like parameter.";
            throw GException::invalid_value(G_GET_MODEL_PARAMETER, msg);
        }
    }

    // ... otherwise find the appropriate model parameter
    else {

        // If a node function was specified then a valid parameter name is
        // needed
        if (m_skymodel->spectral()->type() == "NodeFunction") {
            std::string msg = "Spectral model of source \""+m_srcname+"\" is "
                              "a \"NodeFunction\" but no parameter name was "
                              "specified. Please use the \"parname\" parameter "
                              "to specify the parameter that should be used "
                              "for the upper limit computation.";
            throw GException::invalid_value(G_GET_MODEL_PARAMETER, msg);
        }

        // Find appropriate model parameter and store its pointer in the
        // m_model_par member
        for (const char** par = pars; *par != NULL; ++par) {
            if (m_skymodel->spectral()->has_par(*par)) {
                m_model_par = &(m_skymodel->spectral()->operator[](*par));
                break;
            }
        }

    } // endelse: searched for parameter

    // If no model parameter was found then throw an exception
    if (m_model_par == NULL) {
        std::string msg = "Require one of the following spectral "
                          "parameters for upper limit computation: ";
        for (const char** par = pars; *par != NULL; ++par) {
            if (par != pars) {
                msg.append(", ");
            }
            msg.append("\""+std::string(*par)+"\"");
        }
        msg.append(". The specified source \""+m_srcname+"\" does not "
                   "have such a parameter.");
        throw GException::invalid_value(G_GET_MODEL_PARAMETER, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Determine parameter brackets
 *
 * @param[out] parmin Minimum parameter value.
 * @param[out] parmax Maximum parameter value.
 *
 * @exception GException::invalid_value
 *            Exceeded maximum number of iterations in the bracket search.
 *
 * Determines the brackets for the parameter of interest.
 ***************************************************************************/
void ctulimit::get_parameter_brackets(double& parmin, double& parmax)
{
    // Write header into logger
    log_header1(TERSE, "Find parameter brackets");

    // Extract current value and error
    m_best_value = m_model_par->factor_value();
    m_best_error = m_model_par->factor_error();

    // If parameter error is zero then take parameter value as error
    if (m_best_error == 0.0) {
        m_best_error = m_best_value;
    }

    // Compute parameter bracketing
    parmin = m_best_value - m_sigma_min * m_best_error;
    parmax = m_best_value + m_sigma_max * m_best_error;
    if (m_model_par->has_min() && m_model_par->factor_min() > parmin) {
        parmin = m_model_par->factor_min();
    }
    if (m_model_par->has_max() && m_model_par->factor_max() < parmax) {
        parmax = m_model_par->factor_max();
    }

    // Initialise iteration counter
    int iter = 0;

    // Make sure that upper limit is bracketed
    while (true) {

        // Throw exception if maximum iterations are reached
        if (iter > m_max_iter) {
            std::string msg = "The maximum number of "+gammalib::str(m_max_iter)+
                              " has been reached. You may consider to increase"
                              " the \"max_iter\" parameter and re-run ctulimit.";
            throw GException::invalid_value(G_GET_PARAMETER_BRACKETS, msg);
        }

        // Evaluate logL at upper boundary
        double logL     = evaluate(*m_model_par, parmax);
        double eval_mid = logL - (m_best_logL + m_dlogL);

        // Write current parameter range into logger
        log_value(NORMAL, "Parameter range",
                  "["+gammalib::str(parmin)+", "+gammalib::str(parmax)+"] "
                  "logL("+gammalib::str(parmax)+")="+gammalib::str(logL));

        // If the logL increase is not enough, increase maximum parameter
        // value, otherwise break
        if (eval_mid < 0.0) {
            parmax *= 10.0;
        }
        else {
            break;
        }

        // Increase number of iterations
        iter++;

    } // endwhile

    // Return
    return;
}


/***********************************************************************//**
 * @brief Calculate upper limit using a bisection method
 *
 * @param[in] parmin Minimum parameter value
 * @param[in] parmax Maximum parameter value
 *
 * @exception GException::invalid_value
 *            Exceeded maximum number of iterations in the upper limit
 *            search.
 *
 * Calculates the upper limit using a bisection method. The bisection is
 * stopped if the log-likelihood value is within the tolerance of the target
 * value. The bisection is also stopped putting a warning in the log file if
 * the parameter interval becomes smaller than a tolerance of 1.0e-6. The
 * method stops with an exception if the maximum number of iterations is
 * exhausted.
 ***************************************************************************/
void ctulimit::ulimit_bisection(const double& parmin, const double& parmax)
{
    // Copy values to working values
    double wrk_min = parmin;
    double wrk_max = parmax;

    // Initialise iteration counter and restart flag
    int  iter    = 0;
    bool restart = false;

    // Loop until breaking condition is reached
    while (true) {

        // Throw exception if maximum iterations are reached
        if (iter > m_max_iter) {
            std::string msg = "The maximum number of "+gammalib::str(m_max_iter)+
                              " has been reached. You may consider to increase"
                              " the \"max_iter\" parameter and re-run ctulimit.";
            throw GException::invalid_value(G_UL_BISECTION, msg);
        }

        // Compute center of boundary
        double mid = (wrk_min + wrk_max) / 2.0;

        // Calculate function value
        double logL     = evaluate(*m_model_par, mid);
        double eval_mid = logL - (m_best_logL + m_dlogL);

        // Log information
        log_value(NORMAL, "Iteration "+gammalib::str(iter),
                  "["+gammalib::str(wrk_min)+", "+gammalib::str(wrk_max)+"] "
                  "logL("+gammalib::str(mid)+")="+gammalib::str(logL)+" "
                  "dlogL="+gammalib::str(eval_mid));

        // Check for convergence inside tolerance
        if (std::abs(eval_mid) < m_tol) {
            break;
        }

        // Assess status
        int status = 0;

        // Check whether mid-point is close to best parameter value
        double eps = (mid != 0.0) ? (mid - m_best_value)/mid : mid - m_best_value;
        if (std::abs(eps) < 1.0e-6) {
            std::string msg =
            " *** WARNING: Parameter range mid-point "+gammalib::str(mid)+" is "
            "close to best factor\n"
            "              value "+gammalib::str(m_best_value)+" yet the "
            "log-likelihood value "+gammalib::str(logL)+"\n"
            "              at the mid-point differs significantly from the "
            "log-likelihood value\n"
            "              "+gammalib::str(m_best_logL)+" for the best factor "
            "value.";
            log_string(TERSE, msg);
            status = 1;
        }

        // Check for vanishing parameter range
        eps = (mid != 0.0) ? (wrk_min - wrk_max)/mid : wrk_min - wrk_max;
        if (std::abs(eps) < 1.0e-6) {
            std::string msg =
            " *** WARNING: Parameter range ["+gammalib::str(wrk_min)+", "+
            gammalib::str(wrk_max)+"] has reduced\n"
            "              to a narrow interval without reaching convergence! "
            "The best\n"
            "              log-likelihood value is probably not an absolute "
            "but only a local\n"
            "              minimum.";
            log_string(TERSE, msg);
            status = 2;
        }

        // Restart bisection if a restart is required
        if (!restart && status != 0) {

            // Adopt current logL as best log-likelihood
            restart     = true;
            iter        = 0;
            wrk_min     = parmin;
            wrk_max     = parmax;
            m_best_logL = logL;

            // Signal restart
            std::string msg =
            " Adopt "+gammalib::str(logL)+" as best log-likelihood and restart "
            "bisection";
            log_string(TERSE, msg);

            // Restart
            continue;

        } // endif: restart bisection if required

        // Change boundaries for further iteration
        if (eval_mid > 0.0) {
            wrk_max = mid;     // logL too large, new range = [wrk_min, mid]
        }
        else {
            wrk_min = mid;     // logL too small, new range = [mid, wrk_max]
        }

        // Increment counter
        iter++;

    } // endwhile

    // Return
    return;

}


/***********************************************************************//**
 * @brief Compute upper limit intensity and fluxes
 *
 * Compute upper limit intensity and fluxes, including a correct computation
 * for the diffuse map cube models.
 ***************************************************************************/
void ctulimit::compute_ulimit(void)
{
    // Initialise upper limit intensity and fluxes
    m_diff_ulimit  = 0.0;
    m_flux_ulimit  = 0.0;
    m_eflux_ulimit = 0.0;

    // Get reference energy for differential upper limit
    GEnergy eref = GEnergy(m_eref, "TeV");

    // Create energy range for flux limits
    GEnergy emin = GEnergy(m_emin, "TeV");
    GEnergy emax = GEnergy(m_emax, "TeV");

    // Set pointer to spectral model and initialise spectral nodes model
    // pointer
    GModelSpectral*      spectral = m_skymodel->spectral();
    GModelSpectralNodes* nodes    = NULL;

    // If the spectral model is a diffuse cube then create a node function
    // spectral model that is the product of the diffuse cube node function
    // and the spectral model evaluated at the energies of the node function
    GModelSpatialDiffuseCube* cube =
        dynamic_cast<GModelSpatialDiffuseCube*>(m_skymodel->spatial());
    if (cube != NULL) {


        // Set MC cone to the entire sky. This method call is needed to
        // set-up the cube spectrum
        cube->mc_cone(GSkyRegionCircle(0.0, 0.0, 180.0));

        // Allocate node function to replace the spectral component
        nodes = new GModelSpectralNodes(cube->spectrum());
        for (int i = 0; i < nodes->nodes(); ++i) {
            GEnergy energy    = nodes->energy(i);
            double  intensity = nodes->intensity(i);
            double  value     = spectral->eval(energy);
            nodes->intensity(i, value*intensity);
        }

        // Set the spectral model pointer to the node function. If there are
        // no nodes the spectral pointer is set to NULL since in this case the
        // spectral->flux method will throw an exception
        if (nodes->nodes() > 0) {
            spectral = nodes;
        }
        else {
            spectral = NULL;
        }

    } // endif: spatial model was a diffuse cube

    // Compute upper limit intensity and fluxes
    if (spectral != NULL) {
        m_diff_ulimit  = spectral->eval(eref);
        m_flux_ulimit  = spectral->flux(emin, emax);
        m_eflux_ulimit = spectral->eflux(emin, emax);
    }

    // Free spectral nodes if they have been allocated
    if (nodes != NULL) {
        delete nodes;
    }

    // Return
    return;
}
