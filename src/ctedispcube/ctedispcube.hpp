/***************************************************************************
 *          ctedispcube - Energy dispersion cube generation tool           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Maria Haupt                                      *
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
 * @brief Energy dispersion cube generation tool definition
 * @author Chia-Chun Lu
 */

#ifndef CTEDISPCUBE_HPP
#define CTEDISPCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTEDISPCUBE_NAME    "ctedispcube"
#define CTEDISPCUBE_VERSION "1.1.0"


/***********************************************************************//**
 * @class ctedispcube
 *
 * @brief Energy dispersion cube generation tool
 ***************************************************************************/
class ctedispcube : public ctobservation {

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
    void                 clear(void);
    void                 run(void);
    void                 save(void);
    const GCTACubeEdisp& edispcube(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctedispcube& app);
    void free_members(void);
    void get_parameters(void);
    void init_cube(void);

    // User parameters
    GFilename m_outcube;       //!< Output exposure cube file
    bool      m_addbounds;     //!< Add energies at boundaries?
    GChatter  m_chatter;       //!< Chattiness

    // Protected members
    GCTACubeEdisp m_edispcube; //!< Energy dispersion cube
};


/***********************************************************************//**
 * @brief Return energy dispersion cube
 *
 * @return Energy dispersion cube
 ***************************************************************************/
inline
const GCTACubeEdisp& ctedispcube::edispcube(void) const
{
    return (m_edispcube);
}

#endif /* CTEDISPCUBE_HPP */
