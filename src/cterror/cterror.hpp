/***************************************************************************
 *                 cterror - Parameter error calculation tool              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2016 by Florent Forest                              *
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
#include "ctlikelihood.hpp"

/* __Definitions _________________________________________________________ */
#define CTERROR_NAME    "cterror"
#define CTERROR_VERSION "1.2.0"


/***********************************************************************//**
 * @class cterror
 *
 * @brief Parameter error calculation tool
 ***************************************************************************/
class cterror : public ctlikelihood {

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
    void clear(void);
    void run(void);
    void save(void);

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const cterror& app);
    void   free_members(void);
    void   get_parameters(void);
    double error_bisection(const double& min, const double& max);

    // User parameters
    std::string   m_srcname;      //!< Name of source
    GFilename     m_outmodel;     //!< Output model XML file
    double        m_confidence;   //!< Confidence level
    double        m_tol;          //!< Tolerance for limit determination
    int           m_max_iter;     //!< Maximum number of iterations
    bool          m_apply_edisp;  //!< Apply energy dispersion?
    GChatter      m_chatter;      //!< Chattiness

    // Protected members
    double        m_value;        //!< Parameter value 
    double        m_dlogL;        //!< Likelihood difference for upper limit computation
    GModelPar*    m_model_par;    //!< Pointer to model parameter
    double        m_best_logL;    //!< Best fit log likelihood of given model
};

#endif /* CTERROR_HPP */
