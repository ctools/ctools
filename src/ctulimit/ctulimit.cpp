/***************************************************************************
 *                   ctulimit - Upper limit calculation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Michael Mayer                                    *
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
#define G_GET_PARAMETERS                         "ctulimit::get_parameters()"
#define G_GET_MODEL_PARAMETER               "ctulimit::get_model_parameter()"
#define G_UL_BISECTION             "ctulimit::ul_bisection(double&, double&)"
#define G_EVALUATE                              "ctulimit::evaluate(double&)"

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
ctulimit::ctulimit(void) : ctool(CTULIMIT_NAME, CTULIMIT_VERSION)
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
ctulimit::ctulimit(const GObservations& obs) :
          ctool(CTULIMIT_NAME, CTULIMIT_VERSION)
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
ctulimit::ctulimit(int argc, char *argv[]) :
          ctool(CTULIMIT_NAME, CTULIMIT_VERSION, argc, argv)
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
ctulimit::ctulimit(const ctulimit& app) : ctool(app)
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
 ***************************************************************************/
ctulimit& ctulimit::operator=(const ctulimit& app)
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
void ctulimit::clear(void)
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
 * @brief Computes the upper limit
 *
 * This method calculates the upper limit depending on the given algorithm
 ***************************************************************************/
void ctulimit::run(void)
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

    // Save original models
    GModels models_orig = m_obs.models();

    // Save original log-likelihood. If the value is zero it has never been
    // computed hence we compute it now. 
    m_best_logL = m_obs.logL();
    if (m_best_logL == 0.0) {

        // Write header
        if (logTerse()) {
            log << std::endl;
            log.header1("Compute best-fit likelihood");
        }

        // Reoptimize if likelihood was not given before
        GOptimizerLM* opt = new GOptimizerLM();
        m_obs.optimize(*opt);
        m_obs.errors(*opt);
        m_best_logL = m_obs.logL();

        // Write optimised model into logger
        if (logTerse()) {
            log << m_obs.models() << std::endl;
        }

    } // endif: likelihood was zero

    // Compute parameter bracketing
    double value  = m_model_par->value();
    double error  = m_model_par->error();
    double parmin = value + m_sigma_min * error;
    double parmax = value + m_sigma_max * error;

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Compute upper limit");
        log << gammalib::parformat("Model name");
        log << m_skymodel->name() << std::endl;
        log << gammalib::parformat("Parameter name");
        log << m_model_par->name() << std::endl;
        log << gammalib::parformat("Initial parameter range");
        log << "[";
        log << parmin;
        log << ", ";
        log << parmax;
        log << "]" << std::endl;
    }

    // Compute upper limit
    ulimit_bisection(parmin, parmax);

    // Write final parameter range
    if (logTerse()) {
        log << gammalib::parformat("Final parameter");
        log << m_model_par->value() << std::endl;
    }

    // Get reference energy for differential upper limit
    GEnergy eref = GEnergy(m_eref, "TeV");

    // Create energy range for flux limits
    GEnergy emin = GEnergy(m_emin, "TeV");
    GEnergy emax = GEnergy(m_emax, "TeV");

    // Compute upper limit intensity and fluxes
    m_diff_ulimit  = m_skymodel->spectral()->eval(eref, GTime());
    m_flux_ulimit  = m_skymodel->spectral()->flux(emin, emax);
    m_eflux_ulimit = m_skymodel->spectral()->eflux(emin, emax);

    // Write results to logfile
    if (logTerse()) {
        log << std::endl;
        log.header1("Upper limit results");
        log << gammalib::parformat("Differential flux limit");
        log << m_diff_ulimit;
        log << " ph/cm2/s/MeV at ";
        log << m_eref << " TeV";
        log << std::endl;
        log << gammalib::parformat("Integral flux limit");
        log << m_flux_ulimit;
        log << " ph/cm2/s within [";
        log << m_emin << "-" << m_emax;
        log << "] TeV";
        log << std::endl;
        log << gammalib::parformat("Energy flux limit");
        log << m_eflux_ulimit;
        log << " erg/cm2/s within [";
        log << m_emin << "-" << m_emax;
        log << "] TeV";
        log << std::endl;
    }

    // Recover original models
    m_obs.models(models_orig);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save maps
 *
 * This method saves the upper limit to an ascii file
 ***************************************************************************/
void ctulimit::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save upper limit");
    }

    // Get output filename
    m_outfile = (*this)["outfile"].filename();

    // Create CSV table with 3 columns
    GCsv table(1, 3);

    table.real(0, 0, m_diff_ulimit);
    table.real(0, 1, m_flux_ulimit);
    table.real(0, 2, m_eflux_ulimit);

    // Save CSV table
    table.save(m_outfile, " ", clobber());

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
    m_outfile.clear();
    m_srcname.clear();
    m_sigma_min    = 0.0;
    m_sigma_max    = 0.0;
    m_eref         = 0.0;
    m_emin         = 0.0;
    m_emax         = 0.0;
    m_tol          = 1.0e-6;
    m_max_iter     = 50;

    // Initialise protected members
    m_obs.clear();
    m_models.clear();
    m_dlogL        = 0.0;
    m_skymodel     = NULL;
    m_model_par    = NULL;
    m_best_logL    = 0.0;
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
    m_outfile      = app.m_outfile;
    m_srcname      = app.m_srcname;
    m_sigma_min    = app.m_sigma_min;
    m_sigma_max    = app.m_sigma_max;
    m_eref         = app.m_eref;
    m_emin         = app.m_emin;
    m_emax         = app.m_emax;
    m_tol          = app.m_tol;
    m_max_iter     = app.m_max_iter;


    // Copy protected members
    m_obs          = app.m_obs;
    m_models       = app.m_models;
    m_dlogL        = app.m_dlogL;
    m_best_logL    = app.m_best_logL;
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
 * @exception GException::invalid_value
 *            Test source not found.
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user.
 ***************************************************************************/
void ctulimit::get_parameters(void)
{
    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {

        // Throw exception if no input observation file is given
        require_inobs(G_GET_PARAMETERS);

        // Build observation container
        m_obs = get_observations();

    } // endif: there was no observation in the container


    // If there is are no models associated with the observations then
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

    // Get relevant model and parameter for upper limit computation
    get_model_parameter();

    // Get confidence level
    double confidence = (*this)["confidence"].real();
    if (confidence != 0.95) {
        std::string msg = "Confidence level different from 95% requested."
                          " Currently only 95% is supported.";
        throw GException::invalid_value(G_GET_PARAMETERS, msg);
    }
    else {
        // Set Likelihood difference for 95% CL.
        // See Minuit Handbook
        m_dlogL = 3.84 / 2.0;
    }

    // Read starting boundaries for bisection
    m_sigma_min = (*this)["sigma_min"].real();
    m_sigma_max = (*this)["sigma_max"].real();

    // Read energy values
    m_eref = (*this)["eref"].real();
    m_emin = (*this)["emin"].real();
    m_emax = (*this)["emax"].real();

    // Read precision
    m_tol      = (*this)["tol"].real();
    m_max_iter = (*this)["max_iter"].integer();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outfile = (*this)["outfile"].filename();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * @exception GException::invalid_value
 *            Did not find a valid model parameter
 *
 * Extracts a pointer to the sky model (m_skymodel) and a pointer to the
 * relevant model parameter (m_model_par) from the model container.
 ***************************************************************************/
void ctulimit::get_model_parameter(void)
{
    // Get relevant model and parameter for upper limit computation.
    GModels& models = const_cast<GModels&>(m_obs.models());
    m_skymodel      = dynamic_cast<GModelSky*>(models[m_srcname]);
    if (m_skymodel == NULL) {
        std::string msg = "Source \""+m_srcname+"\" is not a sky model. "
                          "Please specify the name of a sky model for "
                          "upper limit computation.";
        throw GException::invalid_value(G_GET_MODEL_PARAMETER, msg);
    }
    if (m_skymodel->spectral()->has_par("Normalization")) {
        m_model_par = &(m_skymodel->spectral()->operator[]("Normalization"));
    }
    else if (m_skymodel->spectral()->has_par("Prefactor")) {
        m_model_par = &(m_skymodel->spectral()->operator[]("Prefactor"));
    }
    else if (m_skymodel->spectral()->has_par("Integral")) {
        m_model_par = &(m_skymodel->spectral()->operator[]("Integral"));
    }
    else {
        std::string msg = "Require spectral parameter \"Normalization\", "
                          "\"Prefactor\" or \"Integral\" for upper limit "
                          "computation. The specified source \""+m_srcname+
                          "\" does not have such a parameter.";
        throw GException::invalid_value(G_GET_MODEL_PARAMETER, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Performs upper limit computation by using a bisection method
 *
 * @param[in] min Minimum parameter value
 * @param[in] max Maximum parameter value
 *
 * This method calculates the upper limit using a bisection method.
 ***************************************************************************/
void ctulimit::ulimit_bisection(const double& min, const double& max)
{
    // Copy values to working values
    double wrk_min = min;
    double wrk_max = max;

    // Initialise iteration counter
    int iter = 0;

    // Loop until breaking condition is reached
    while (true) {

        // Log information
        if (logExplicit()) {
            log << gammalib::parformat("Iteration "+gammalib::str(iter));
            log << "[";
            log << wrk_min;
            log << ", ";
            log << wrk_max;
            log << "]" << std::endl;
        }

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
        double eval_mid = evaluate(mid);

        // Check for convergence inside tolerance
        if (std::abs(eval_mid) < m_tol) {
            break;
        }

        // Change boundaries for further iteration
        if (eval_mid > 0.0) {
            wrk_max = mid;
        }
        else if (eval_mid < 0.0) {
            wrk_min = mid;
        }

        // Increment counter
        iter++;

    } // endwhile

    // Return
    return;

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
double ctulimit::evaluate(const double& value)
{
    // Initialise log-likelihood value
    double logL = 0.0;

    // Check if given parameter is within boundaries
    if (value > m_model_par->min() && value < m_model_par->max()) {

        // Change parameter factor
        m_model_par->value(value);

        // Evaluate likelihood for new model container
        m_obs.eval();

        // Retrieve likelihood
        logL = m_obs.logL();

    } // endif: value was inside allowed range

    // ... otherwise signal that the parameter went outside the boundaries
    else {
        std::string msg = "Value of parameter \""+m_model_par->name()+"\" "
                          "outside of validity range requested. To omit "
                          "this error please enlarge the  parameter range "
                          "in the model XML file.";
        throw GException::invalid_value(G_EVALUATE, msg);
    }

    // Compute function value
    double logL_difference = logL - m_best_logL - m_dlogL;

    // Return
    return logL_difference;
}
