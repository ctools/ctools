/***************************************************************************
 *           ctprob - Computes probability for a given model               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Leonardo Di Venere                          *
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
 * @file ctprob.hpp
 * @brief Computes probability for a given model
 * @author Leonardo Di Venere
 */

#ifndef CTPROB_HPP
#define CTPROB_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTPROB_NAME    "ctprob"
#define CTPROB_VERSION "0.0.1"


/***********************************************************************//**
 * @class ctprob
 *
 * @brief Data selection tool
 ***************************************************************************/
class ctprob : public ctobservation {

public:
    // Constructors and destructors
    ctprob(void);
    explicit ctprob(const GObservations& obs);
    ctprob(int argc, char *argv[]);
    ctprob(const ctprob& app);
    virtual ~ctprob(void);

    // Operators
    ctprob& operator=(const ctprob& app);

    // Methods
    void clear(void);
    void run(void);
    void save(void);
    void publish(const std::string& name = "");

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const ctprob& app);
    void        free_members(void);
    void        get_parameters(void);
    void        get_obs(void);
    void        evaluate_probability(GCTAObservation*   obs);
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
    std::string m_outobs;     //!< Output event list or XML file
    std::string m_prefix;     //!< Prefix for multiple event lists
    bool        m_apply_edisp;  //!< Apply energy dispersion?
    bool        m_publish;      //!< Publish model cube?
    GChatter    m_chatter;    //!< Chattiness

    // Protected members
    std::vector<std::string> m_infiles;       //!< Input event filenames
    std::vector<std::string> m_evtname;       //!< Event extension names
    std::vector<std::string> m_gtiname;       //!< GTI extension names
};


#endif /* CTPROB_HPP */
