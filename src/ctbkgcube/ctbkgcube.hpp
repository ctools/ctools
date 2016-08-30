/***************************************************************************
 *               ctbkgcube - Background cube generation tool               *
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
 * @file ctbkgcube.hpp
 * @brief Background cube generation tool definition
 * @author Chia-Chun Lu
 */

#ifndef CTBKGCUBE_HPP
#define CTBKGCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTBKGCUBE_NAME    "ctbkgcube"
#define CTBKGCUBE_VERSION "1.1.0"


/***********************************************************************//**
 * @class ctbkgcube
 *
 * @brief Background cube generation tool
 ***************************************************************************/
class ctbkgcube : public ctobservation {

public:
    // Constructors and destructors
    ctbkgcube(void);
    explicit ctbkgcube(const GObservations& obs);
    ctbkgcube(int argc, char *argv[]);
    ctbkgcube(const ctbkgcube& app);
    virtual ~ctbkgcube(void);

    // Operators
    ctbkgcube& operator=(const ctbkgcube& app);

    // Methods
    void                      clear(void);
    void                      run(void);
    void                      save(void);
    void                      publish(const std::string& name = "");
    const GCTACubeBackground& bkgcube(void) const;
    const GModels&            models(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctbkgcube& app);
    void free_members(void);
    void get_parameters(void);
    void init_cube(void);

    // User parameters
    GFilename          m_outcube;     //!< Filename of output cube
    GFilename          m_outmodel;    //!< Filename of output XML model
    bool               m_addbounds;   //!< Add energies at boundaries?
    bool               m_publish;     //!< Publish background cube?
    GChatter           m_chatter;     //!< Chattiness

    // Protected members
    GCTACubeBackground m_background;  //!< Background cube
    GModels            m_bkgmdl;      //!< CTA background models
    GModels            m_outmdl;      //!< Output models
};


/***********************************************************************//**
 * @brief Return background response cube containing background rate
 *
 * @return Background response cube.
 ***************************************************************************/
inline
const GCTACubeBackground& ctbkgcube::bkgcube(void) const
{
    // Return background cube
    return (m_background);
}



/***********************************************************************//**
 * @brief Return background model container
 *
 * @return Background model container.
 ***************************************************************************/
inline
const GModels& ctbkgcube::models(void) const
{
    // Return background model container
    return (m_outmdl);
}

#endif /* CTBKGCUBE_HPP */
