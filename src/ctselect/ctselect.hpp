/***************************************************************************
 *                      ctselect - Data selection tool                     *
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
 * @file ctselect.hpp
 * @brief Data selection tool definition
 * @author Juergen Knoedlseder
 */

#ifndef CTSELECT_HPP
#define CTSELECT_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctool.hpp"

/* __Definitions _________________________________________________________ */
#define CTSELECT_NAME    "ctselect"
#define CTSELECT_VERSION "00-07-00"


/***********************************************************************//**
 * @class ctselect
 *
 * @brief Data selection tool
 ***************************************************************************/
class ctselect : public ctool {

public:
    // Constructors and destructors
    ctselect(void);
    explicit ctselect(const GObservations& obs);
    ctselect(int argc, char *argv[]);
    ctselect(const ctselect& app);
    virtual ~ctselect(void);

    // Operators
    ctselect& operator=(const ctselect& app);

    // Methods
    void                  clear(void);
    void                  run(void);
    void                  save(void);
    const GObservations&  obs(void) const;
    const GTimeReference& time_reference(void) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const ctselect& app);
    void        free_members(void);
    void        get_parameters(void);
    void        select_events(GCTAObservation* obs, const std::string& filename);
    GEbounds    set_ebounds(GCTAObservation* obs,
                            const GEbounds& ebounds) const;
    std::string check_infile(const std::string& filename) const;
    std::string set_outfile_name(const std::string& filename) const;
    void        save_fits(void);
    void        save_xml(void);
    void        save_event_list(const GCTAObservation* obs,
                                const std::string&     infile,
                                const std::string&     outfile) const;

    // User parameters
    std::string m_outobs;    //!< Output event list or XML file
    std::string m_prefix;     //!< Prefix for multiple event lists
    bool        m_usepnt;     //!< Use pointing instead of RA/DEC parameters
    double      m_ra;         //!< RA of ROI centre
    double      m_dec;        //!< DEC of ROI centre
    double      m_rad;        //!< ROI radius
    double      m_tmin;       //!< Start time
    double      m_tmax;       //!< Stop time
    double      m_emin;       //!< Lower energy
    double      m_emax;       //!< Upper energy
    std::string m_expr;       //!< Selection expression
    std::string m_usethres;   //!< Energy threshold type

    // Protected members
    GObservations            m_obs;        //!< Observations container
    std::vector<std::string> m_infiles;    //!< Input event filenames
    GTime                    m_timemin;    //!< Earliest time
    GTime                    m_timemax;    //!< Latest time

};


/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GObservations& ctselect::obs(void) const
{
    return m_obs;
}


/***********************************************************************//**
 * @brief Return CTA time reference
 *
 * @return CTA time reference
 ***************************************************************************/
inline
const GTimeReference& ctselect::time_reference(void) const
{
    return m_cta_ref;
}

#endif /* CTSELECT_HPP */
