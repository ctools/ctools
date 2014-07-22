/***************************************************************************
 *                     ctpsfcube - CTA PSF cube tool                       *
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
 * @file ctpsfcube.hpp
 * @brief CTA PSF cube tool definition
 * @author Chia-Chun Lu
 */

#ifndef CTPSFCUBE_HPP
#define CTPSFCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTPSFCUBE_NAME    "ctpsfcube"
#define CTPSFCUBE_VERSION "00-01-00"


/***********************************************************************//**
 * @class ctpsfcube
 *
 * @brief CTA PSF cube tool
 ***************************************************************************/
class ctpsfcube : public GApplication  {

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
    void                 clear(void);
    void                 execute(void);
    void                 run(void);
    void                 save(void);
    const GObservations& obs(void) const;
    const GCTAMeanPsf&   psfcube(void) const;
    void                 get_parameters(void);
    void                 get_obs(void);
    void                 set_response(void);
    void                 get_ebounds(void);
    void                 set_from_cntmap(const std::string& filename);

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const ctpsfcube& app);
    void           free_members(void);

    // User parameters
    std::string   m_outfile;     //!< Output PSF cube file
    bool          m_apply_edisp; //!< Apply energy dispersion?

    // Protected members
    bool          m_read_ahead;  //!< Read ahead parameters
    GObservations m_obs;         //!< Observation container
    GCTAMeanPsf   m_psfcube;     //!< PSF cube
    GEbounds      m_ebounds;     //!< Energy boundaries
};


/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Observation container
 ***************************************************************************/
inline
const GObservations& ctpsfcube::obs(void) const
{
    return (m_obs);
}


/***********************************************************************//**
 * @brief Return PSF cube
 *
 * @return PSF cube
 ***************************************************************************/
inline
const GCTAMeanPsf& ctpsfcube::psfcube(void) const
{
    return (m_psfcube);
}

#endif /* CTPSFCUBE_HPP */
