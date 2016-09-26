/***************************************************************************
 *                    cttsmap - TS map calculation tool                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Michael Mayer                               *
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
 * @brief TS map calculation tool interface definition
 * @author Michael Mayer
 */

#ifndef CTTSMAP_HPP
#define CTTSMAP_HPP

/* __ Includes ___________________________________________________________ */
#include "ctlikelihood.hpp"

/* __Definitions _________________________________________________________ */
#define CTTSMAP_NAME    "cttsmap"
#define CTTSMAP_VERSION "1.2.0"


/***********************************************************************//**
 * @class cttsmap
 *
 * @brief TS map calculation tool
 *
 * This class computes a set of maps from any kind of observations:
 *
 *      - TS map
 *      - Flux map
 *      - Map of spectral index
 *
 * The class operates on predefined observation containers, an individual
 * event list or an observation definition XML file.
 *
 * During the computation a putative point-like source is moved along a
 * grid of coordinates. The best fit results (TS, flux, index) are stored in
 * maps which are saved in the output FITS files.
 ***************************************************************************/
class cttsmap : public ctlikelihood {

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
    void           clear(void);
    void           run(void);
    void           save(void);
    void           publish(const std::string& name = "");
    const GSkyMap& tsmap(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const cttsmap& app);
    void free_members(void);
    void get_parameters(void);
    void init_maps(const GSkyMap& map);

    // User parameters
    std::string m_srcname;     //!< Name of source which is moved around
    GFilename   m_outmap;      //!< Output TS map file name
    bool        m_apply_edisp; //!< Apply energy dispersion?
    bool        m_publish;     //!< Publish TS map?
    bool        m_errors;      //!< Compute and store parameter errors?

    // Parameters to control speed and job splitting
    int         m_binmin;      //!< Map bin number from which computation should start
    int         m_binmax;      //!< Map bin number where map computation should end
    double      m_logL0;       //!< Likelihood value of null hypothesis

    // Protected members
    GSkyMap                  m_tsmap;       //!< TS map
    GSkyMap                  m_statusmap;   //!< Map of optimizer fit status
    std::vector<std::string> m_mapnames;    //!< Names of free parameters
    std::vector<GSkyMap>     m_maps;        //!< Sky maps for each free parameter
    GModel*                  m_testsource;  //!< Pointer to test source for TS computation
};


/***********************************************************************//**
 * @brief Return TS skymap
 *
 * @return Reference to TS skymap
 ***************************************************************************/
inline
const GSkyMap& cttsmap::tsmap(void) const
{
    return m_tsmap;
}

#endif /* CTTSMAP_HPP */
