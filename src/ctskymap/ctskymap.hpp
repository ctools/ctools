/***************************************************************************
 *                       ctskymap - Sky mapping tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
#define CTSKYMAP_VERSION "1.2.0"


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
    void           publish(const std::string& name = "");
    const GSkyMap& map(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctskymap& app);
    void free_members(void);
    void get_parameters(void);
    void map_events(GCTAObservation* obs);

    // User parameters
    GFilename     m_outmap;     //!< Output file name
    double        m_emin;       //!< Minimum energy (TeV)
    double        m_emax;       //!< Maximum energy (TeV)
    bool          m_publish;    //!< Publish sky map?
    GChatter      m_chatter;    //!< Chattiness

    // Protected members
    GObservations m_obs;        //!< Observation container
    GSkyMap       m_skymap;     //!< Sky map
};


/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GSkyMap& ctskymap::map(void) const
{
    return m_skymap;
}

#endif /* CTSKYMAP_HPP */
