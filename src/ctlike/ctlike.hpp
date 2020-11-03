/***************************************************************************
 *                ctlike - Maximum likelihood fitting tool                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
    void clear(void);
    void run(void);
    void save(void);

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const ctlike& app);
    void   free_members(void);
    void   get_parameters(void);
    void   optimize_lm(void);
    double reoptimize_lm(void);
    GXml   xml_result(void) const;

    // User parameters
    GFilename m_outmodel;        //!< Source model output XML file name
    GFilename m_outcovmat;       //!< Covariance matrix output file name
    bool      m_refit;           //!< Refitting
    bool      m_apply_edisp;     //!< Apply energy dispersion?
    bool      m_fix_spat_for_ts; //!< Fix spatial parameters for TS computation?
    GChatter  m_chatter;         //!< Chattiness

    // Members
    double    m_logL;            //!< Maximum log likelihood
    double    m_nobs;            //!< Number of observed events
    double    m_npred;           //!< Number of predicted events
};

#endif /* CTLIKE_HPP */
