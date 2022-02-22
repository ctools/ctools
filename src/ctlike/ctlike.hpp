/***************************************************************************
 *                ctlike - Maximum likelihood fitting tool                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2022 by Juergen Knoedlseder                         *
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
 * @file ctlike.hpp
 * @brief Maximum likelihood fitting tool definition
 * @author Juergen Knoedlseder
 */

#ifndef CTLIKE_HPP
#define CTLIKE_HPP

/* __ Includes ___________________________________________________________ */
#include "ctlikelihood.hpp"

/* __Definitions _________________________________________________________ */
#define CTLIKE_NAME "ctlike"


/***********************************************************************//**
 * @class ctlike
 *
 * @brief Maximum likelihood fitting tool
 ***************************************************************************/
class ctlike : public ctlikelihood {

public:
    // Constructors and destructors
    ctlike(void);
    explicit ctlike(const GObservations& obs);
    ctlike(int argc, char *argv[]);
    ctlike(const ctlike& app);
    virtual ~ctlike(void);

    // Operators
    ctlike& operator=(const ctlike& app);

    // Methods
    void          clear(void);
    void          process(void);
    void          save(void);
    const int&    iter(void) const;
    const double& logL(void) const;
    const double& nobs(void) const;
    const double& npred(void) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const ctlike& app);
    void   free_members(void);
    void   get_parameters(void);
    void   optimize_lm(void);
    double reoptimize_lm(void);
    GXml   xml_result(void) const;
    bool   refit(const GOptimizer* opt);

    // User parameters
    GFilename m_outmodel;        //!< Source model output XML file name
    GFilename m_outcovmat;       //!< Covariance matrix output file name
    int       m_max_iter;        //!< Maximum number of iterations
    bool      m_refit;           //!< Refitting?
    bool      m_refit_if_failed; //!< Refitting in case of failure?
    bool      m_apply_edisp;     //!< Apply energy dispersion?
    bool      m_fix_spat_for_ts; //!< Fix spatial parameters for TS computation?
    GChatter  m_chatter;         //!< Chattiness

    // Members
    int       m_iter;            //!< Number of iterations
    double    m_logL;            //!< Maximum log likelihood
    double    m_nobs;            //!< Number of observed events
    double    m_npred;           //!< Number of predicted events
};


/***********************************************************************//**
 * @brief Return number of maximum likelihood iterations
 *
 * @return Number of maximum likelihood iterations.
 ***************************************************************************/
inline
const int& ctlike::iter(void) const
{
    return (m_iter);
}


/***********************************************************************//**
 * @brief Return maximum likelihood value
 *
 * @return Maximum likelihood value.
 ***************************************************************************/
inline
const double& ctlike::logL(void) const
{
    return (m_logL);
}


/***********************************************************************//**
 * @brief Return number of observed events
 *
 * @return Number of observed events.
 ***************************************************************************/
inline
const double& ctlike::nobs(void) const
{
    return (m_nobs);
}


/***********************************************************************//**
 * @brief Return number of predicted events
 *
 * @return Number of predicted events.
 ***************************************************************************/
inline
const double& ctlike::npred(void) const
{
    return (m_npred);
}

#endif /* CTLIKE_HPP */
