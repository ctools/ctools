/***************************************************************************
 *              ctlikelihood - Base class for likelihood tools             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2021 by Juergen Knoedlseder                         *
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
 * @file ctlikelihood.cpp
 * @brief Likelihood tool base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <typeinfo>
#include "ctlikelihood.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPT                                "ctlikelihood::opt(GOptimizer*)"
#define G_EVALUATE              "ctlikelihood::evaluate(GModelPar&, double&)"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Name constructor
 *
 * @param[in] name Likelihood tool name.
 * @param[in] version Likelihood tool version.
 *
 * Constructs a likelihood tool from the @p name and @p version. See the
 * equivalent ctool constructor for details.
 ***************************************************************************/
ctlikelihood::ctlikelihood(const std::string& name,
                           const std::string& version) :
                           ctobservation(name, version)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Application parameters constructor
 *
 * @param[in] name Observation tool name.
 * @param[in] version Observation tool version.
 * @param[in] pars Application parameters.
 *
 * Constructs a likelihood tool from the @p name, @p version and the
 * application parameters @p pars. See the equivalent ctool constructor
 * for details.
 ***************************************************************************/
ctlikelihood::ctlikelihood(const std::string&      name,
                           const std::string&      version,
                           const GApplicationPars& pars) :
                           ctobservation(name, version, pars)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Command line constructor
 *
 * @param[in] name Likelihood tool name.
 * @param[in] version Likelihood tool version.
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 *
 * Constructs a likelihood tool from the @p name, @p version and command
 * line arguments. See the equivalent ctool constructor for details.
 ***************************************************************************/
ctlikelihood::ctlikelihood(const std::string& name,
                           const std::string& version,
                           int   argc,
                           char *argv[]) :
                           ctobservation(name, version, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observations constructor
 *
 * @param[in] name Likelihood tool name.
 * @param[in] version Likelihood tool version.
 * @param[in] obs Observation container.
 *
 * Constructs a likelihood tool from the @p name, @p version and an
 * observation container.
 ***************************************************************************/
ctlikelihood::ctlikelihood(const std::string&   name,
                           const std::string&   version,
                           const GObservations& obs) :
                           ctobservation(name, version, obs)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Likelihood tool.
 *
 * Constructs an instance of a likelihood tool by copying information from
 * another likelihood tool.
 ***************************************************************************/
ctlikelihood::ctlikelihood(const ctlikelihood& app) : ctobservation(app)
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
 * Destructs the likelihood tool.
 ***************************************************************************/
ctlikelihood::~ctlikelihood(void)
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
 * @param[in] app Likelihood tool.
 * @return Likelihood tool.
 *
 * Assigns a likelihood tool.
 ***************************************************************************/
ctlikelihood& ctlikelihood::operator=(const ctlikelihood& app)
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
 * @brief Set optimizer
 *
 * @param[in] opt Pointer to optimizer.
 *
 * @exception GException::invalid_argument
 *            Specified optimizer is not of type GOptimizerLM.
 *
 * Sets the optimizer. So far only optimizers of type GOptimizerLM are
 * supported. If a different optimizer is specified the method throws an
 * exception.
 ***************************************************************************/
void ctlikelihood::opt(const GOptimizer* opt)
{
    // Throw an exception if the optimizer is not a LM optimizer
    const GOptimizerLM* lm = dynamic_cast<const GOptimizerLM*>(opt);
    if (lm == NULL) {
        std::string cls = std::string(typeid(opt).name());
        std::string msg = "Method requires \"GOptimizerLM\" optimizer yet the "
                          "provided optimizer is of type \""+cls+"\". Please "
                          "specify an instance of \"GOptimizerLM\" as "
                          "argument.";
        throw GException::invalid_argument(G_OPT, msg);
    }

    // Set optimizer
    m_opt = *lm;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                   Protected methods exposed in Python                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Evaluates the log-likelihood function
 *
 * @param[in] par Model parameter
 * @param[in] value Model parameter factor value
 * @return Log-likelihood function
 *
 * Evaluates the log-likelihood function at a given @p value.
 ***************************************************************************/
double ctlikelihood::evaluate(GModelPar& par, const double& value)
{
    // Initialise log-likelihood value
    double logL = 0.0;

    // Throw an exception if the parameter is below the minimum boundary
    if (par.has_min() && value < par.factor_min()) {
        std::string msg = "Parameter \""+par.name()+"\" value"+
                          gammalib::str(value)+" is below its minimum boundary "+
                          gammalib::str(par.factor_min())+". To omit this "
                          "error please lower the minimum parameter boundary.";
        throw GException::invalid_value(G_EVALUATE, msg);
    }

    // Throw an exception if the parameter is above the maximum boundary
    if (par.has_max() && value > par.factor_max()) {
        std::string msg = "Parameter \""+par.name()+"\" value"+
                          gammalib::str(value)+" is above its maximum boundary "+
                          gammalib::str(par.factor_max())+". To omit this "
                          "error please raise the maximum parameter boundary.";
        throw GException::invalid_value(G_EVALUATE, msg);
    }

    // Change parameter factor
    par.factor_value(value);

    // Fix parameter
    par.fix();

    // Re-optimize log-likelihood
    m_obs.optimize(m_opt);

    // Free parameter
    par.free();

    // Retrieve log-likelihood
    logL = m_obs.logL();

    // Return log-likelihood
    return logL;
}


/*==========================================================================
 =                                                                         =
 =                            Protected methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void ctlikelihood::init_members(void)
{
    // Initialise members
    m_opt.clear();

    // Initialise optimizer characteristics
    m_opt.max_iter(50);     // Maximum number of iterations
    m_opt.max_stalls(10);   // Maximum number of stalls
    m_opt.eps(0.005);       // Accuracy of maximum likelihood value

    // Set optimizer characteristics from optional user parameters
    if (has_par("like_accuracy")) {
        m_opt.eps((*this)["like_accuracy"].real());
    }
    if (has_par("max_iter")) {
        m_opt.max_iter((*this)["max_iter"].integer());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Likelihood tool.
 ***************************************************************************/
void ctlikelihood::copy_members(const ctlikelihood& app)
{
    // Copy members
    m_opt = app.m_opt;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctlikelihood::free_members(void)
{
    // Return
    return;
}
