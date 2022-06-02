/***************************************************************************
 *                      ctselect - Data selection tool                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2022 by Juergen Knoedlseder                         *
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
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTSELECT_NAME "ctselect"


/***********************************************************************//**
 * @class ctselect
 *
 * @brief Data selection tool
 ***************************************************************************/
class ctselect : public ctobservation {

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
    void clear(void);
    void process(void);
    void save(void);
    void publish(const std::string& name = "");

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const ctselect& app);
    void        free_members(void);
    void        get_parameters(void);
    void        select_events(GCTAObservation* obs,
                              const std::string& filename,
                              const std::string& evtname,
                              const std::string& gtiname);
    GEbounds    set_ebounds(GCTAObservation* obs,
                            const GEbounds& ebounds) const;
    std::string check_infile(const std::string& filename,
                             const std::string& evtname) const;
    void        save_fits(void);
    void        save_xml(void);

    // User parameters
    std::string m_outobs;     //!< Output event list or XML file
    double      m_emin;       //!< Lower energy
    double      m_emax;       //!< Upper energy
    std::string m_expr;       //!< Selection expression
    std::string m_usethres;   //!< Energy threshold type
    GChatter    m_chatter;    //!< Chattiness
    bool        m_forcesel;   //!< Enforce RoI selection

    // Protected members
    std::vector<std::string> m_infiles;       //!< Input event filenames
    std::vector<std::string> m_evtname;       //!< Event extension names
    std::vector<std::string> m_gtiname;       //!< GTI extension names
    GPhases                  m_phases;        //!< Phase intervals
    bool                     m_select_energy; //!< Perform energy selection
    bool                     m_select_phase;  //!< Perform phase selection
};

#endif /* CTSELECT_HPP */
