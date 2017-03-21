/***************************************************************************
 *                       ctskymap - Sky mapping tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2017 by Juergen Knoedlseder                         *
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
#include <string>
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTSKYMAP_NAME    "ctskymap"
#define CTSKYMAP_VERSION "1.3.0"


/***********************************************************************//**
 * @class ctskymap
 *
 * @brief Sky mapping tool
 *
 * This class creates a sky map from a CTA event list.
 ***************************************************************************/
class ctskymap : public ctobservation {

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
    const GSkyMap& skymap(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctskymap& app);
    void free_members(void);
    void get_parameters(void);
    void map_events(GCTAObservation* obs);
    void map_background(GCTAObservation* obs);
    void map_background_irf(GCTAObservation* obs);

    // Background integration kernel
    class irf_kern : public GFunction {
    public:
        irf_kern(const GCTABackground* bgd,
                 const GCTAInstDir*    dir) :
                 m_bgd(bgd),
                 m_dir(dir) { }
        double eval(const double& lnE);
    protected:
        const GCTABackground* m_bgd;   //!< Pointer to background
        const GCTAInstDir*    m_dir;   //!< Pointer to instrument direction
    };

    // User parameters
    GFilename     m_outmap;      //!< Output file name
    double        m_emin;        //!< Minimum energy (TeV)
    double        m_emax;        //!< Maximum energy (TeV)
    std::string   m_bkgsubtract; //!< Background subtraction method
    bool          m_publish;     //!< Publish sky map?
    GChatter      m_chatter;     //!< Chattiness

    // Protected members
    GSkyMap       m_skymap;     //!< Sky map
    GSkyMap       m_bkgmap;     //!< Background map
    GSkyMap       m_sigmap;     //!< Significance map
};


/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GSkyMap& ctskymap::skymap(void) const
{
    return m_skymap;
}

#endif /* CTSKYMAP_HPP */
