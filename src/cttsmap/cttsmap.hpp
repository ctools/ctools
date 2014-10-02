/***************************************************************************
 *                      cttsmap - TS map tool                      *
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
 * @file cttsmap.hpp
 * @brief TS map calculation tool
 * @author Michael Mayer
 */

#ifndef CTTSMAP_HPP
#define CTTSMAP_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTTSMAP_NAME    "cttsmap"
#define CTTSMAP_VERSION "00-01-00"


/***********************************************************************//**
 * @class cttsmap
 *
 * @brief TS map calculation interface
 *
 * This class computes a set of maps from any kind of observations:
 *  - TS map
 *  - flux map
 *  - map of spectral index
 *  The class operates on predefined observation containers, an individual
 *   event list or an observation definition XML file.
 *
 *  During the computation a putative point-like source is moved along a
 *  grid of coordinates. The best fit results (TS, flux, index) are stored in maps
 *  which are saved in the output FITS files.
 *
 ***************************************************************************/
class cttsmap : public GApplication  {
public:
    // Constructors and destructors
	cttsmap(void);
    explicit cttsmap(const GObservations& obs);
    cttsmap(int argc, char *argv[]);
    cttsmap(const cttsmap& app);
    virtual ~cttsmap(void);

    // Operators
    cttsmap& operator=(const cttsmap& app);

    // Methods
    void                 clear(void);
    void                 execute(void);
    void                 run(void);
    void                 save(void);
    const GObservations& obs(void) const;
    const GSkymap& tsmap(void) const;
    const GSkymap& fluxmap(void) const;
    const GSkymap& indexmap(void) const;
    const GSkymap& statusmap(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const cttsmap& app);
    void free_members(void);
    void get_parameters(void);
    void init_maps(void);
    void init_testsource(void);

    // User parameters
    std::string              m_evfile;     //!< Input event list or XML file
    std::string              m_outfile;    //!< Output counts map or XML file
    std::string              m_proj;       //!< WCS projection
    std::string              m_coordsys;   //!< Coordinate system
    double                   m_xref;       //!< Longitude reference coordinate
    double                   m_yref;       //!< Latitude reference coordinate
    double                   m_binsz;      //!< Pixel size
    int                      m_nxpix;      //!< Number of pixels in longitude
    int                      m_nypix;      //!< Number of pixels in latitude

    // Parameters to control speed and job splitting
    int                      m_binmin;  //!< Map bin number from which computation should start
    int                      m_binmax;  //!< Map bin number where map computation should end
    double              m_logL0; //!< Likelihood value of null hypothesis

    // Protected members
    GObservations            m_obs;        //!< Observation container
    bool                     m_read_ahead; //!< Read ahead parameters
    GSkymap                  m_tsmap;       //!< ts map
    GSkymap                  m_fluxmap;    //!< flux map
    GSkymap                  m_indexmap; //!< spectral index map
    GSkymap                  m_statusmap; //!< map of computed bins
    GModelSky               m_testsource; //!< test source for TS computation
    double               m_source_integral; //!< Integral flux of test source
    double			   m_source_integral_min; //!< minimum boundary of integral flux
    double			   m_source_integral_max; //!< maximum boundary of integral flux
    double               m_source_index; //!< assumed spectral index of test source
    double               m_source_index_min; //!< minimum boundary of spectral index of test source
    double               m_source_index_max; //!< maximum boundary of spectral index of test source
    GEnergy             m_source_emin; //!< Minimum energy for source integral
    GEnergy             m_source_emax; //!< Maximum energy for source integral
    bool                   m_free_index; //!< flag if spectral index of test source is left free or fixed

};

/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GObservations& cttsmap::obs(void) const
{
    return m_obs;
}

/***********************************************************************//**
 * @brief Return TS skymap
 *
 * @return Reference to TS skymap
 ***************************************************************************/
inline
const GSkymap& cttsmap::tsmap(void) const
{
    return m_tsmap;
}

/***********************************************************************//**
 * @brief Return flux skymap
 *
 * @return Reference to flux skymap
 ***************************************************************************/
inline
const GSkymap& cttsmap::fluxmap(void) const
{
    return m_fluxmap;
}

/***********************************************************************//**
 * @brief Return spectral index skymap
 *
 * @return Reference to spectral index skymap
 ***************************************************************************/
inline
const GSkymap& cttsmap::indexmap(void) const
{
    return m_indexmap;
}

/***********************************************************************//**
 * @brief Return status skymap
 *
 * @return Reference to status skymap
 ***************************************************************************/
inline
const GSkymap& cttsmap::statusmap(void) const
{
    return m_statusmap;
}


#endif /* CTTSMAP_HPP */
