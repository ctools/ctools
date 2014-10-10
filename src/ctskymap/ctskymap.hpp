/***************************************************************************
 *                       ctskymap - Sky mapping tool                       *
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
 * @file ctskymap.hpp
 * @brief Sky mapping tool definition
 * @author Juergen Knoedlseder
 */

#ifndef CTSKYMAP_HPP
#define CTSKYMAP_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctool.hpp"

/* __Definitions _________________________________________________________ */
#define CTSKYMAP_NAME    "ctskymap"
#define CTSKYMAP_VERSION "00-02-00"


/***********************************************************************//**
 * @class ctskymap
 *
 * @brief Sky mapping tool
 *
 * This class creates a sky map from a CTA event list.
 ***************************************************************************/
class ctskymap : public ctool {

public:
    // Constructors and destructors
    ctskymap(void);
    explicit ctskymap(const GObservations& obs);
    ctskymap(int argc, char *argv[]);
    ctskymap(const ctskymap& app);
    virtual ~ctskymap(void);

    // Operators
    ctskymap& operator=(const ctskymap& app);

    // Methods
    void           clear(void);
    void           run(void);
    void           save(void);
    const GSkymap& map(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctskymap& app);
    void free_members(void);
    void get_parameters(void);
    void init_map(GCTAObservation* obs);
    void map_events(GCTAObservation* obs);

    // User parameters
    std::string   m_evfile;     //!< Input event list
    std::string   m_outfile;    //!< Output counts map
    double        m_emin;       //!< Minimum energy (TeV)
    double        m_emax;       //!< Maximum energy (TeV)
    std::string   m_proj;       //!< WCS projection
    std::string   m_coordsys;   //!< Coordinate system
    double        m_xref;       //!< Longitude reference coordinate
    double        m_yref;       //!< Latitude reference coordinate
    double        m_binsz;      //!< Pixel size
    int           m_nxpix;      //!< Number of pixels in longitude
    int           m_nypix;      //!< Number of pixels in latitude

    // Protected members
    GObservations m_obs;        //!< Observation container
    GSkymap       m_skymap;     //!< Sky map
};


/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GSkymap& ctskymap::map(void) const
{
    return m_skymap;
}

#endif /* CTSKYMAP_HPP */
