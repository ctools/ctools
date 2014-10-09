/***************************************************************************
 *                ctmodel - CTA model cube generation tool                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @file ctmodel.hpp
 * @brief CTA model cube tool definition
 * @author Juergen Knoedlseder
 */

#ifndef CTMODEL_HPP
#define CTMODEL_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTMODEL_NAME    "ctmodel"
#define CTMODEL_VERSION "00-02-00"


/***********************************************************************//**
 * @class ctmodel
 *
 * @brief CTA counts model tool interface defintion
 *
 * This class creates counts model map(s). The definition of the counts model
 * can be taken from a predefined observation container, from the counts maps
 * found in an observation definition XML file, from an individual counts
 * map, or from the parameter definition.
 *
 * Results are stored in an observation container that can be written to disk
 * in form of FITS files (model maps) and an updated observation definition
 * XML file.
 ***************************************************************************/
class ctmodel : public GApplication  {

public:
    // Constructors and destructors
    ctmodel(void);
    explicit ctmodel(const GObservations& obs);
    ctmodel(int argc, char *argv[]);
    ctmodel(const ctmodel& app);
    virtual ~ctmodel(void);

    // Operators
    ctmodel& operator=(const ctmodel& app);

    // Methods
    void                 clear(void);
    void                 execute(void);
    void                 run(void);
    void                 save(void);
    const GObservations& obs(void) const;
    const GCTAEventCube& cube(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctmodel& app);
    void free_members(void);
    void get_parameters(void);
    void setup_obs(void);
    void init_cube(void);
    void fill_cube(const GCTAObservation* obs);
    void get_ebounds(void);
    
    // User parameters
    std::string m_infile;      //!< Input counts cube
    std::string m_obsfile;     //!< Event list, counts map or observation definition XML file
    std::string m_outfile;     //!< Output model cube
    std::string m_caldb;       //!< Calibration database
    std::string m_irf;         //!< Instrument response functions
    std::string m_srcmdl;      //!< Source model
    double      m_ra;          //!< RA of pointing direction
    double      m_dec;         //!< DEC of pointing direction
    double      m_deadc;       //!< Deadtime correction factor
    double      m_tmin;        //!< Start time
    double      m_tmax;        //!< Stop time
    std::string m_ebinalg;     //!< Algorithm for energy binning
    std::string m_ebinfile;    //!< FITS-file containing energy binning
    double      m_emin;        //!< Lower energy
    double      m_emax;        //!< Upper energy
    int         m_enumbins;    //!< Number of energy bins
    std::string m_proj;        //!< WCS projection
    std::string m_coordsys;    //!< Coordinate system
    double      m_xref;        //!< Longitude reference coordinate
    double      m_yref;        //!< Latitude reference coordinate
    double      m_binsz;       //!< Pixel size
    int         m_nxpix;       //!< Number of pixels in longitude
    int         m_nypix;       //!< Number of pixels in latitude
    bool        m_apply_edisp; //!< Apply energy dispersion?

    // Protected members
    GObservations m_obs;        //!< Observation container
    GCTAEventCube m_cube;       //!< Model cube
    GEbounds      m_ebounds;    //!< Energy boundaries
    GGti          m_gti;        //!< Model cube GTIs
    bool          m_read_ahead; //!< Read ahead parameters
};


/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GObservations& ctmodel::obs(void) const
{
    return m_obs;
}


/***********************************************************************//**
 * @brief Return model cube
 *
 * @return Reference to model cube
 ***************************************************************************/
inline
const GCTAEventCube& ctmodel::cube(void) const
{
    return m_cube;
}

#endif /* CTMODEL_HPP */
