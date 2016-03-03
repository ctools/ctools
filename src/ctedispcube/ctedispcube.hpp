/***************************************************************************
 *                  ctedispcube - EDISP cube generation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Chia-Chun Lu                                *
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
 * @file ctedispcube.hpp
 * @brief EDISP cube generation tool definition
 * @author Chia-Chun Lu
 */

#ifndef CTEDISPCUBE_HPP
#define CTEDISPCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctool.hpp"

/* __Definitions _________________________________________________________ */
#define CTEDISPCUBE_NAME    "ctedispcube"
#define CTEDISPCUBE_VERSION "1.1.0"


/***********************************************************************//**
 * @class ctedispcube
 *
 * @brief EDISP cube generation tool
 ***************************************************************************/
class ctedispcube : public ctool {

public:
    // Constructors and destructors
    ctedispcube(void);
    explicit ctedispcube(const GObservations& obs);
    ctedispcube(int argc, char *argv[]);
    ctedispcube(const ctedispcube& app);
    virtual ~ctedispcube(void);

    // Operators
    ctedispcube& operator=(const ctedispcube& app);

    // Methods
    void               clear(void);
    void               run(void);
    void               save(void);
    const GCTACubeEdisp& edispcube(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctedispcube& app);
    void free_members(void);
    void get_parameters(void);

    // User parameters
    GFilename     m_outcube;     //!< Output EDISP cube file name
    bool          m_apply_edisp; //!< Apply energy dispersion?

    // Protected members
    GObservations m_obs;         //!< Observation container
    GCTACubeEdisp   m_edispcube;     //!< EDISP cube
};


/***********************************************************************//**
 * @brief Return EDISP cube
 *
 * @return EDISP cube
 ***************************************************************************/
inline
const GCTACubeEdisp& ctedispcube::edispcube(void) const
{
    return (m_edispcube);
}

#endif /* CTEDISPCUBE_HPP */
