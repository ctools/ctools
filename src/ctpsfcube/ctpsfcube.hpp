/***************************************************************************
 *                  ctpsfcube - PSF cube generation tool                   *
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
 * @file ctpsfcube.hpp
 * @brief PSF cube generation tool definition
 * @author Chia-Chun Lu
 */

#ifndef CTPSFCUBE_HPP
#define CTPSFCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTPSFCUBE_NAME    "ctpsfcube"
#define CTPSFCUBE_VERSION "1.1.0"


/***********************************************************************//**
 * @class ctpsfcube
 *
 * @brief PSF cube generation tool
 ***************************************************************************/
class ctpsfcube : public ctobservation {

public:
    // Constructors and destructors
    ctpsfcube(void);
    explicit ctpsfcube(const GObservations& obs);
    ctpsfcube(int argc, char *argv[]);
    ctpsfcube(const ctpsfcube& app);
    virtual ~ctpsfcube(void);

    // Operators
    ctpsfcube& operator=(const ctpsfcube& app);

    // Methods
    void               clear(void);
    void               run(void);
    void               save(void);
    const GCTACubePsf& psfcube(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctpsfcube& app);
    void free_members(void);
    void get_parameters(void);
    void init_cube(void);

    // User parameters
    GFilename m_outcube;         //!< Output PSF cube file name
    bool      m_addbounds;       //!< Add energies at boundaries?
    GChatter  m_chatter;         //!< Chattiness

    // Protected members
    GCTACubePsf   m_psfcube;     //!< PSF cube
};


/***********************************************************************//**
 * @brief Return PSF cube
 *
 * @return PSF cube
 ***************************************************************************/
inline
const GCTACubePsf& ctpsfcube::psfcube(void) const
{
    return (m_psfcube);
}

#endif /* CTPSFCUBE_HPP */
