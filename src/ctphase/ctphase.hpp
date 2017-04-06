/***************************************************************************
 *          ctphase - Append phase information to CTA events file          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Leonardo Di Venere                               *
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
 * @file ctphase.hpp
 * @brief Append phase information to CTA events file
 * @author Leonardo Di Venere
 */

#ifndef CTPHASE_HPP
#define CTPHASE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTPHASE_NAME    "ctphase"
#define CTPHASE_VERSION "0.0.1"


/***********************************************************************//**
 * @class ctphase
 *
 * @brief Append phase columm to observation files
 ***************************************************************************/
class ctphase : public ctobservation {

public:
    // Constructors and destructors
    ctphase(void);
    explicit ctphase(const GObservations& obs);
    ctphase(int argc, char *argv[]);
    ctphase(const ctphase& app);
    virtual ~ctphase(void);

    // Operators
    ctphase& operator=(const ctphase& app);
    
    // Methods
    void clear(void);
    void run(void);
    void save(void);
    void publish(const std::string& name = "");

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const ctphase& app);
    void        free_members(void);
    void        get_parameters(void);
    void        read_phase_info_from_xml(void);
    void        phase_events(GCTAObservation* obs,
                             const std::string& filename,
                             const std::string& evtname,
                             const std::string& gtiname);
    std::string check_infile(const std::string& filename,
                             const std::string& evtname) const;
    std::string set_outfile_name(const std::string& filename) const;
    std::string get_gtiname(const std::string& filename,
                            const std::string& evtname) const;
    void        save_fits(void);
    void        save_xml(void);
    void        save_event_list(const GCTAObservation* obs,
                                const std::string&     infile,
                                const std::string&     evtname,
                                const std::string&     gtiname,
                                const std::string&     outfile) const;

    // User parameters
    std::string m_outobs;               //!< Output event list or XML file
    std::string m_prefix;               //!< Prefix for multiple event lists
    GModelTemporalPhaseCurve m_phase;   //!< Phase information object
    GChatter    m_chatter;              //!< Chattiness

    // Protected members
    std::vector<std::string> m_infiles;       //!< Input event filenames
    std::vector<std::string> m_evtname;       //!< Event extension names
    std::vector<std::string> m_gtiname;       //!< GTI extension names
};


#endif /* CTPHASE_HPP */
