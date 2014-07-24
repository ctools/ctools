/***************************************************************************
 *                   ctbkgcube - CTA background cube tool                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Chia-Chun Lu                                     *
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
 * @brief CTA background cube tool definition
 * @author Chia-Chun Lu
 */

#ifndef CTBKGCUBE_HPP
#define CTBKGCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTBKGCUBE_NAME    "ctbkgcube"
#define CTBKGCUBE_VERSION "00-01-00"


/***********************************************************************//**
 * @class ctbkgcube
 *
 * @brief CTA background cube tool
 ***************************************************************************/
class ctbkgcube : public GApplication  {

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
    void                 clear(void);
    void                 execute(void);
    void                 run(void);
    void                 save(void);
    const GCTAEventCube& cube(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctbkgcube& app);
    void free_members(void);
    void get_parameters(void);
    void get_obs(void);
    void get_ebounds(void);
    void set_from_cntmap(const std::string& filename);
    void fill_cube(GCTAObservation* obs);

    // Protected members
    bool          m_read_ahead;  //!< Read ahead parameters
    GObservations m_obs;         //!< Observation container
    GCTAEventCube m_bkgcube;     //!< Background cube
    GModels       m_bkgmdl;      //!< CTA background models
    GEbounds      m_ebounds;     //!< Energy boundaries
    double        m_livetime;    //!< Total livetime for normalization
};


/***********************************************************************//**
 * @brief Return CTA event cube containing background rate
 *
 * @return CTA event cube.
 ***************************************************************************/
inline
const GCTAEventCube& ctbkgcube::cube(void) const
{
    // Return background cube
    return (m_bkgcube);
}

#endif /* CTBKGCUBE_HPP */
