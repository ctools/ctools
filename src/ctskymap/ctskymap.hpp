/***************************************************************************
 *                       ctskymap - Sky mapping tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2019 by Juergen Knoedlseder                         *
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
    const GFits&   fits(void) const;
    void           exclusion_map(const GSkyRegionMap& exclusion_map);
    GSkyRegionMap  exclusion_map(void) const;

protected:
    // Protected methods
    void     init_members(void);
    void     copy_members(const ctskymap& app);
    void     free_members(void);
    void     get_parameters(void);
    void     setup_maps(void);
    void     setup_exclusion_map(void);
    void     setup_exclusion_map_fits(const GFilename& filename);
    void     setup_exclusion_map_region(const GFilename& filename);
    void     adjust_exclusion_map(void);
    void     fill_maps(void);
    void     fill_maps_counts(GCTAObservation* obs);
    void     fill_maps_acceptance(GCTAObservation* obs);
    void     compute_maps(void);
    void     compute_maps_ring_fft(void);
    void     compute_maps_ring_direct(void);
    void     compute_ring_values(const int&     ipixel,
                                 const GSkyMap& counts,
                                 const GSkyMap& background,
                                 double&        non,
                                 double&        noff,
                                 double&        alpha);
    void     construct_fits(void);
    void     ring_bounding_box(const int& ipixel, int& ix1, int& ix2,
                                                  int& iy1, int& iy2) const;
    GSkyMap  ring_convolve(const GSkyMap& map, const double& rmin,
                                               const double& rmax) const;
    GNdarray ring_kernel(const double& rmin, const double& rmax) const;
    double   sigma_li_ma(const double& n_on,
                         const double& n_off,
                         const double& alpha) const;
    void     write_map(GFits&             fits,
                       const GSkyMap&     map,
                       const std::string& extname) const;
    void     write_hdu_keywords(GFitsHDU* hdu) const;

    // User parameters
    GFilename     m_outmap;      //!< Output file name
    GFilename     m_inexclusion; //!< Exclusion map file name
    double        m_emin;        //!< Minimum energy (TeV)
    double        m_emax;        //!< Maximum energy (TeV)
    std::string   m_bkgsubtract; //!< Background subtraction method
    double        m_roiradius;   //!< Region of interest radius for RING bkg.
    double        m_inradius;    //!< Inner ring radius for RING background
    double        m_outradius;   //!< Outer ring radius for RING background
    int           m_iterations;  //!< Number of iterations for RING background
    double        m_threshold;   //!< Threshold for RING background
    bool          m_usefft;      //!< Use FFT for RING background
    bool          m_publish;     //!< Publish sky map?
    GChatter      m_chatter;     //!< Chattiness

    // Protected members
    bool          m_has_inmap;   //!< Has valid input map
    GSkyMap       m_skymap;      //!< Sky map
    GSkyMap       m_bkgmap;      //!< Background map
    GSkyMap       m_sigmap;      //!< Significance map
    GSkyMap       m_counts;      //!< Counts map
    GSkyMap       m_acceptance;  //!< Acceptance map
    GSkyMap       m_exclmap;     //!< Exclusion map for RING background
    GFits         m_fits;        //!< Output GFits object

    // Caching variables to prevent multiple computations of the same thing
    std::vector<double>  m_solidangle;    //!< Cached pixel solid angles
    std::vector<GSkyDir> m_dirs;          //!< Cached pixel directions
    double               m_cos_roiradius; //!< Cosine of RoI radius
    double               m_cos_inradius;  //!< Cosine of inner ring radius
    double               m_cos_outradius; //!< Cosine of outer ring radius
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


/***********************************************************************//**
 * @brief Return fits container
 *
 * @return Reference to fits container
 ***************************************************************************/
inline
const GFits& ctskymap::fits(void) const
{
    return m_fits;
}


/***********************************************************************//**
 * @brief Set exclusion region map
 *
 * @param[in] exclusion_map Exclusion region map.
 ***************************************************************************/
inline
void ctskymap::exclusion_map(const GSkyRegionMap& exclusion_map)
{
    // Assign exclusion map
    m_exclmap = exclusion_map.map();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return exclusion region map
 *
 * @return Exclusion region map
 ***************************************************************************/
inline
GSkyRegionMap ctskymap::exclusion_map(void) const
{
    GSkyRegionMap exclusion_map(m_exclmap);
    return exclusion_map;
}

#endif /* CTSKYMAP_HPP */
