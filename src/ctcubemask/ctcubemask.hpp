/***************************************************************************
 *                      ctcubemask - Cube filter tool                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Chia-Chun Lu                                *
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
 * @brief Cube filter tool definition
 * @author Chia-Chun Lu
 */

#ifndef CTCUBEMASK_HPP
#define CTCUBEMASK_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTCUBEMASK_NAME    "ctcubemask"
#define CTCUBEMASK_VERSION "1.2.0"


/***********************************************************************//**
 * @class ctcubemask
 *
 * @brief Cube filter tool
 ***************************************************************************/
class ctcubemask : public ctobservation {

public:
    // Constructors and destructors
    ctcubemask(void);
    explicit ctcubemask(const GObservations& obs);
    ctcubemask(int argc, char *argv[]);
    ctcubemask(const ctcubemask& app);
    virtual ~ctcubemask(void);

    // Operators
    ctcubemask& operator=(const ctcubemask& app);

    // Methods
    void clear(void);
    void run(void);
    void save(void);
    void publish(const std::string& name = "");

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const ctcubemask& app);
    void        free_members(void);
    void        get_parameters(void);
    void        apply_mask(GCTAObservation* obs);
    std::string region_string(const GSkyRegion& region) const;
    std::string set_outfile_name(const std::string& filename) const;
    void        save_fits(void);
    void        save_xml(void);

    // User parameters
	GFilename   m_regfile;    //!< ds9 region file
    GFilename   m_outcube;    //!< Output event list or XML file
	std::string m_prefix;     //!< Prefix for multiple counts maps
    bool        m_usepnt;     //!< Use pointing instead of RA/DEC parameters
    GCTARoi     m_roi;        //!< RoI selection
    double      m_emin;       //!< Lower energy
    double      m_emax;       //!< Upper energy
    bool        m_publish;    //!< Publish counts cube?

    // Protected members
    std::vector<std::string> m_infiles;       //!< Input event filenames
    bool                     m_select_energy; //!< Perform energy selection
    bool                     m_select_roi;    //!< Perform ROI selection
};


#endif /* CTCUBEMASK_HPP */
