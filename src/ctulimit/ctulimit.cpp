/***************************************************************************
 *                    ctulimit - upper limit calculation tool                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Michael Mayer                                    *
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
 * @brief upper limit calculation tool interface implementation
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
#define G_GET_PARAMETERS                          "ctulimit::get_parameters()"
#define G_UL_BISECTION        "ctulimit::ul_bisection(const double&, const double&)"
#define G_EVALUATE        "ctulimit::evaluate(const double&)"
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

    GModels models_orig = m_obs.models();


    m_bestloglike = m_obs.logL();
    if (m_bestloglike == 0.0) {

        // Write header
        if (logTerse()) {
            log << std::endl;
            log.header1("Recompute best-fit likelihood");
        }

        // Reoptimize if likelihood was not given before
        GOptimizerLM* opt = new GOptimizerLM();
        m_obs.optimize(*opt);
        m_obs.errors(*opt);
        m_bestloglike = m_obs.logL();

    } //endif: likelihood was 0.0

    // Store optimized models and pointer to skymodel of interest
    m_models = m_obs.models();
    m_skymodel = dynamic_cast<GModelSky*>(m_models[m_srcname])->clone();

    // Get reference energy for differential upper limit
    GEnergy eref = GEnergy(m_eref, "TeV");

    // Create energy range for flux limits
    GEnergy emin = GEnergy(m_emin, "TeV");
    GEnergy emax = GEnergy(m_emax, "TeV");

    // Get value, scale and error of parameter
    double value = (*m_skymodel->spectral())[0].value();
    double error = (*m_skymodel->spectral())[0].error();
    double scale = (*m_skymodel->spectral())[0].scale();

    // Compute starting boundaries
    double parmin = value / scale +  m_sigma_min * error / scale;
    double parmax = value / scale +  m_sigma_max * error / scale;

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Compute upper limit");
        log << "Searching for upper limit between ";
        log << parmin * scale;
        log << " and ";
        log << parmax * scale;
        log << ": "<<std::endl;
    }

    // compute upper limit
    ulimit_bisection(parmin, parmax, scale);

    m_diff_ulimit = m_skymodel->spectral()->eval(eref, GTime());
    m_flux_ulimit = m_skymodel->spectral()->flux(emin, emax);
    m_eflux_ulimit = m_skymodel->spectral()->eflux(emin, emax);

    // Write results to logfile
    if (logTerse()) {
        log << std::endl;
        log.header1("Upper limit computation finished");
        log << "Upper limits are:" << std::endl;
        log << " Differential at ";
        log << m_eref <<" TeV: ";
        log << m_diff_ulimit << " ph/cm2/s/MeV" << std::endl;
        log << " Integral [";
        log << m_emin << " - " << m_emax;
        log << "] TeV: ";
        log << m_flux_ulimit << " ph/cm2/s" <<std::endl;
        log << " Energy flux [";
        log << m_emin << " - " << m_emax;
        log << "] TeV: ";
        log << m_eflux_ulimit << " erg/cm2/s" <<std::endl;
    }

    // Recover original models
    m_obs.models(models_orig);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Evaluates the likelihood
 *
 * @param[in] value Factorised parameter value of the parameter of interest
 * @return likelihood value
 *
 * This method evaluates the likelihood as a function of the parameter of interest
 ***************************************************************************/
double ctulimit::evaluate(const double& value)
{
    // Initialise likelihood value
    double LogL = 0.0;

    // Check if given parameter is within boundaries
    if (value > (*m_skymodel->spectral())[0].factor_min() && value < (*m_skymodel->spectral())[0].factor_max()) {

        // Remove source model
        m_models.remove(m_srcname);

        // Change spectral model value
        (*m_skymodel->spectral())[0].factor_value(value);

        // Re-append model to container
        m_models.append(*m_skymodel);

        // Assign new models to observations
        m_obs.models(m_models);

        // Evaluate likelihood for new model container
        m_obs.eval();

        // Retrieve likelihood
        LogL = m_obs.logL();

    } // endif: value was inside allowed range

    else {
        // throw exception if value out of range
        std::string msg = "Value out of range requested: "
                "To omit this error, you could increase the allowed "
                "parameter range in the model xml file.";
        throw GException::invalid_value(G_EVALUATE, msg);
    }

    // Compute function value
    double difflike = LogL - m_bestloglike - m_dloglike;

    // Return
    return difflike;
}

/***********************************************************************//**
 * @brief Induces upper limit computation by using a bisection method
 *
 * This method calculates the upper limit using a bisection method
 ***************************************************************************/
void ctulimit::ulimit_bisection(const double& min, const double& max, const double& scale)
{
    // copy values to working values
    double wrk_min = min;
    double wrk_max = max;

    // Initialise counter
    int iter=0;

    // Loop until breaking condition is reached
    while (true) {

        // Log information
            if (logExplicit()) {
                log << " Iteration ";
                log << iter;
                log << ": Parameter range reduced to [";
                log << wrk_min * scale;
                log <<"; ";
                log << wrk_max * scale;
                log << "]" << std::endl;
            }

        // Throw exception if maximum iterations are reached
        if( iter > m_max_iter) {
            throw GException::invalid_value(G_UL_BISECTION, "Maximum iterations reached");;
        }

        // compute center of boundary
        double mid = (wrk_min + wrk_max) / 2.0;

        // Calculate function value
        double eval_mid = evaluate(mid);

        // Check for convergence inside tolerance
        if (std::abs(eval_mid) < m_tol) {
            break;
        }

        // change boundaries for further iteration
        if (eval_mid > 0.0) {
            wrk_max = mid;
        }
        else if (eval_mid < 0.0) {
            wrk_min = mid;
        }

        // increment counter
        iter++;

    } // endwhile

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
    // Initialise members
    m_outfile.clear();
    m_srcname.clear();
    m_bestloglike = 0.0;
    m_flux_ulimit = 0.0;
    m_diff_ulimit = 0.0;
    m_eflux_ulimit = 0.0;
    m_dloglike = 0.0;
    m_tol = 1e-6;
    m_max_iter = 50;
    m_eref = 0.0;
    m_emin = 0.0;
    m_emax = 0.0;
    m_sigma_min = 0.0;
    m_sigma_max  = 0.0;
    m_skymodel = NULL;

    // Initialise protected members
    m_obs.clear();
    m_models.clear();

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
    // Copy attributes
    m_outfile  = app.m_outfile;
    m_srcname = app.m_srcname;
    m_diff_ulimit = app.m_diff_ulimit;
    m_flux_ulimit = app.m_flux_ulimit;
    m_eflux_ulimit = app.m_eflux_ulimit;
    m_bestloglike = app.m_bestloglike;
    m_dloglike = app.m_dloglike;
    m_tol = app.m_tol;
    m_max_iter = app.m_max_iter;
    m_eref = app.m_eref;
    m_emin = app.m_emin;
    m_emax = app.m_emax;
    m_sigma_min = app.m_sigma_min;
    m_sigma_max = app.m_sigma_max;

    // Copy protected members
    m_obs        = app.m_obs;
    m_models  = app.m_models;

    // Clone protected members
    m_skymodel = (app.m_skymodel != NULL) ? app.m_skymodel->clone() : NULL;

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
 *            Test source not found or no RA/DEC parameters found for test
 *            source.
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. Most parameters are only required if no observation exists so
 * far in the observation container. In this case, a single CTA observation
 * will be added to the container, using the definition provided in the
 * parameter file.
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

    double CL = (*this)["cl"].real();
    if (CL != 0.95) {
        std::string msg = "Confidence level different from 95% requested"
                          "Currently only 95% CL is possible";
        throw GException::invalid_value(G_GET_PARAMETERS, msg);
    }
    else {
        // Set Likelihood difference for 95% CL.
        // See Minuit Handbook
        m_dloglike = 3.84 / 2.0;
    }

    // Read starting boundaries for bisection
    m_sigma_min = (*this)["sigma_min"].real();
    m_sigma_max = (*this)["sigma_max"].real();

    // Read energy values
    m_eref = (*this)["eref"].real();
    m_emin = (*this)["emin"].real();
    m_emax = (*this)["emax"].real();

    // Read precision
    m_tol = (*this)["tol"].real();
    m_max_iter = (*this)["max_iter"].integer();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outfile = (*this)["outfile"].filename();
    }

    // Return
    return;
}


