/***************************************************************************
 *                   ctlike - CTA maximum likelihood tool                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @brief CTA maximum likelihood tool interface definition
 * @author Juergen Knoedlseder
 */

#ifndef CTLIKE_HPP
#define CTLIKE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTLIKE_NAME    "ctlike"
#define CTLIKE_VERSION "00-05-00"


/***********************************************************************//**
 * @class ctlike
 *
 * @brief CTA maximum likelihood tool interface defintion.
 ***************************************************************************/
class ctlike : public GApplication  {
public:
    // Constructors and destructors
    ctlike(void);
    explicit ctlike(const GObservations& obs);
    ctlike(int argc, char *argv[]);
    ctlike(const ctlike& app);
    virtual ~ctlike(void);

    // Operators
    ctlike& operator= (const ctlike& app);

    // Methods
    void                 clear(void);
    void                 execute(void);
    void                 run(void);
    void                 save(void);
    const GObservations& obs(void) const;
    const GOptimizer*    opt(void) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const ctlike& app);
    void   free_members(void);
    void   get_parameters(void);
    void   optimize_lm(void);
    double reoptimize_lm(void);

    // User parameters
    std::string   m_stat;        //!< Optimisation statistics (poisson/gaussian)
    bool          m_refit;       //!< Refitting
    bool          m_tscalc;      //!< TS computation
    std::string   m_caldb;       //!< Calibration database
    std::string   m_irf;         //!< Instrument response functions
    std::string   m_outmdl;      //!< Source model output XML file
    bool          m_apply_edisp; //!< Apply energy dispersion?

    // Members
    GObservations m_obs;        //!< Observations
    int           m_max_iter;   //!< Maximum number of iterations
    int           m_max_stall;  //!< Maximum number of stalls
    double        m_logL;       //!< Maximum log likelihood
    GOptimizer*   m_opt;        //!< Optimizer
    bool          m_read_ahead; //!< Read ahead parameters
};


/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GObservations& ctlike::obs(void) const
{
    return m_obs;
}


/***********************************************************************//**
 * @brief Return optimizer
 *
 * @return Pointer to optimizer
 ***************************************************************************/
inline
const GOptimizer* ctlike::opt(void) const
{
    return m_opt;
}

#endif /* CTLIKE_HPP */
