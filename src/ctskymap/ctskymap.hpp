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
#define CTSKYMAP_NAME "ctskymap"


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
    void map_background_ring(GCTAObservation* obs);
    void map_significance(void);
    void map_exclusions(const GFilename& filename);
    void map_exclusions_fits(const GFilename& filename);
    void map_exclusions_reg(const GFilename& filename);
    void compute_ring_values(const GSkyMap& counts, 
                             const GSkyMap& sensitivity,
                             const GSkyDir& position,
                             double& non, double& noff, double& alpha);
    void write_hdu_keywords(GFitsHDU* hdu) const;

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
    double        m_roiradius;   //!< Region of interest radius for RING bkg.
    double        m_inradius;    //!< Inner ring radius for RING bkg.
    double        m_outradius;   //!< Outer ring radius for RING bkg.
    bool          m_publish;     //!< Publish sky map?
    GChatter      m_chatter;     //!< Chattiness

    // Protected members
    GSkyMap       m_skymap;      //!< Sky map
    GSkyMap       m_bkgmap;      //!< Background map
    GSkyMap       m_sigmap;      //!< Significance map
    GSkyMap       m_exclmap;     //!< Exclusion map for RING background
    GSkyMap       m_alphamap;    //!< Alpha map values for RING background
    GSkyMap       m_onmap;       //!< On counts map for RING background

    // Caching variables to prevent multiple computations of the same thing
    std::vector<double>  m_solidangle; //!< Cached pixel solid angles
    std::vector<GSkyDir> m_dirs;       //!< Cached pixel directions
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
