/***************************************************************************
 *                    ctulimit - Upper limit calculation tool                    *
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
 * @brief upper limit calculation tool interface definition
 * @author Michael Mayer
 */

#ifndef CTULIMIT_HPP
#define CTULIMIT_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctool.hpp"

/* __Definitions _________________________________________________________ */
#define CTULIMIT_NAME    "ctulimit"
#define CTULIMIT_VERSION "00-01-00"


/***********************************************************************//**
 * @class ctulimit
 *
 * @brief Upper limit calculation tool
 *
 * This class computes and upper limit for a given parameter.
 *
 * The class operates on predefined observation containers, an individual
 * event list or an observation definition XML file.
 *
 * During the computation the likelihood function is inspected to find the best
 * value for the upper limit.
 *
 ***************************************************************************/
class ctulimit : public ctool {

public:
    // Constructors and destructors
    ctulimit(void);
    explicit ctulimit(const GObservations& obs);
    ctulimit(int argc, char *argv[]);
    ctulimit(const ctulimit& app);
    virtual ~ctulimit(void);

    // Operators
    ctulimit& operator=(const ctulimit& app);

    // Methods
    void                 clear(void);
    void                 run(void);
    void                 save(void);
    const GObservations& obs(void) const;
    const double&     diff_ulimit(void) const;
    const double&     flux_ulimit(void) const;
    const double&    eflux_ulimit(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctulimit& app);
    void free_members(void);
    void get_parameters(void);
    void ulimit_bisection(const double& parmin, const double& parmax, const double& scale);
    double evaluate(const double& value);

    // User parameters
    std::string              m_srcname;    //!< Name of source which is moved around
    std::string              m_outfile;    //!< Output ascii file

    GModelSky*            m_skymodel; //!< Pointer to variable spectral model

    int                m_max_iter; //!< maximum number of iterations
    double                m_dloglike; //!< Likelihood difference for upper limit computation
    double                m_bestloglike ; //!< Best fit log likelihood of given model
    double                m_tol; //!< tolerance for limit determination
    double                m_sigma_min; //!< Starting value minimum (multiple fit errors above fit values)
    double                m_sigma_max; //!< Starting value maximum (multiple fit errors above fit values)
    double                m_eref; //!< reference energy for flux limits
    double                m_emin; //!< minimum energy for flux limits
    double                m_emax; //!< maximum energy for flux limits

    double                m_diff_ulimit; //!<Differential upper limit value
    double                m_flux_ulimit; //!< Flux upper limit value
    double                m_eflux_ulimit; //!< Energy flux upper limits

    // Protected members
    GObservations            m_obs;        //!< Observation container
    GModels m_models;

};


/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GObservations& ctulimit::obs(void) const
{
    return m_obs;
}

/***********************************************************************//**
 * @brief Return differential upper limit
 *
 * @return differential upper limit
 ***************************************************************************/
inline
const double& ctulimit::diff_ulimit(void) const
{
    return m_diff_ulimit;
}

/***********************************************************************//**
 * @brief Return flux upper limit
 *
 * @return flux upper limit
 ***************************************************************************/
inline
const double& ctulimit::flux_ulimit(void) const
{
    return m_flux_ulimit;
}

/***********************************************************************//**
 * @brief return energy flux upper limit
 *
 * @return energy flux upper limits
 ***************************************************************************/
inline
const double& ctulimit::eflux_ulimit(void) const
{
    return m_eflux_ulimit;
}

#endif /* CTULIMIT_HPP */
