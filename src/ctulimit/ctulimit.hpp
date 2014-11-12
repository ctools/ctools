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
    const double&     ulimit(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctulimit& app);
    void free_members(void);
    void get_parameters(void);

    // User parameters
    std::string              m_infile;     //!< Input file
    std::string              m_modelfile;//!< Name of model optimised model file
    std::string              m_srcname;    //!< Name of source which is moved around
    std::string              m_parname; //!< Name of parameter upper limit should be computed
    std::string              m_outfile;    //!< Output counts map or XML file
    std::string              m_algorithm; //!< algorithm of upper limit computation

    // Protected members
    GObservations            m_obs;        //!< Observation container

    double                m_ulimit; //!<Upper limit value

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
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const double& ctulimit::ulimit(void) const
{
    return m_ulimit;
}


#endif /* CTULIMIT_HPP */
