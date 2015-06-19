/***************************************************************************
 *                   cterror - Parameter error calculation tool            *
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
 * @file cterror.hpp
 * @brief Parameter error calculation tool interface definition
 * @author Florent Forest
 */

#ifndef CTERROR_HPP
#define CTERROR_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctool.hpp"

/* __Definitions _________________________________________________________ */
#define CTERROR_NAME    "cterror"
#define CTERROR_VERSION "01-00-00"


/***********************************************************************//**
 * @class cterror
 *
 * @brief Parameter error calculation tool
 *
 ***************************************************************************/
class cterror : public ctool {

public:
    // Constructors and destructors
    cterror(void);
    explicit cterror(const GObservations& obs);
    cterror(int argc, char *argv[]);
    cterror(const cterror& app);
    virtual ~cterror(void);

    // Operators
    cterror& operator=(const cterror& app);

    // Methods
    void                 clear(void);
    void                 run(void);
    void                 save(void);
    const GObservations& obs(void) const;
    const GOptimizer*    opt(void) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const cterror& app);
    void   free_members(void);
    void   get_parameters(void);
    double evaluate(const double& value);
    void   error_bisection(const double& min, const double& max);

    // User parameters
    std::string   m_srcname;      //!< Name of source which is moved around
    double        m_confidence;   //!< Confidence level
    double        m_sigma_min;    //!< Starting value minimum (multiple fit errors above fit values)
    double        m_sigma_max;    //!< Starting value maximum (multiple fit errors above fit values)
    double        m_tol;          //!< Tolerance for limit determination
    int           m_max_iter;     //!< Maximum number of iterations
    double        m_value;        //!< Parameter value 
    double        m_error;        //!< Parameter error

    // Protected members
    GObservations m_obs;          //!< Observation container
    double        m_dlogL;        //!< Likelihood difference for upper limit computation
    GModelSky*    m_skymodel;     //!< Pointer to sky model
    GModelPar*    m_model_par;    //!< Pointer to model parameter
    double        m_best_logL;    //!< Best fit log likelihood of given model
    GOptimizerLM* m_opt;

};


/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GObservations& cterror::obs(void) const
{
    return m_obs;
}

/***********************************************************************//**
 * @brief Return optimizer
 *
 * @return Pointer to optimizer
 ***************************************************************************/
inline
const GOptimizer* cterror::opt(void) const
{
    return m_opt;
}

#endif /* CTERROR_HPP */
