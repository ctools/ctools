/***************************************************************************
 *                    ctcubemask - CTA cube filter tool                      *
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
 * @file ctcubemask.hpp
 * @brief CTA mask cube tool definition
 * @author Chia-Chun Lu
 */

#ifndef CTCUBEMASK_HPP
#define CTCUBEMASK_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTCUBEMASK_NAME    "ctcubemask"
#define CTCUBEMASK_VERSION "00-05-00"


/***********************************************************************//**
 * @class ctcubemask
 *
 * @brief CTA mask cube tool interface defintion.
 ***************************************************************************/
class ctcubemask : public GApplication  {
public:
    // Constructors and destructors
    ctcubemask(void);
    explicit ctcubemask(GObservations obs);
    ctcubemask(int argc, char *argv[]);
    ctcubemask(const ctcubemask& app);
    virtual ~ctcubemask(void);

    // Operators
    ctcubemask& operator= (const ctcubemask& app);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           apply_mask(GCTAObservation* obs, const std::string& filename);

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const ctcubemask& app);
    void        free_members(void);
    std::string check_infile(const std::string& filename) const;
    std::string set_outfile_name(const std::string& filename) const;
    void        save_fits(void);
    void        save_xml(void);
    void        save_counts_map(const GCTAObservation* obs,
                                const std::string&     outfile) const;

    // User parameters
    std::string m_infile;     //!< Input event list or XML file
	std::string m_regfile;    //!< ds9 region file
    std::string m_outfile;    //!< Output event list or XML file
	std::string m_prefix;     //!< Prefix for multiple counts maps
    bool        m_usepnt;     //!< Use pointing instead of RA/DEC parameters
    double      m_ra;         //!< RA of ROI centre
    double      m_dec;        //!< DEC of ROI centre
    double      m_rad;        //!< ROI radius
    double      m_emin;       //!< Lower energy
    double      m_emax;       //!< Upper energy

    // Protected members
    GObservations            m_obs;        //!< Observations container
    std::vector<std::string> m_infiles;    //!< Input event filenames
    bool                     m_use_xml;    //!< Use XML file instead of FITS file
    bool                     m_read_ahead; //!< Read ahead parameters
};

#endif /* CTCUBEMASK_HPP */
