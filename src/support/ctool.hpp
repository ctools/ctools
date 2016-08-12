/***************************************************************************
 *                        ctool - ctool base class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Juergen Knoedlseder                         *
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
 * @file ctool.hpp
 * @brief ctool base class definition
 * @author Juergen Knoedlseder
 */

#ifndef CTOOL_HPP
#define CTOOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */


/***********************************************************************//**
 * @class ctool
 *
 * @brief Base class for ctools
 *
 * This class is the generic base class for all ctools. It enforces
 * implementation of the clear(), run() and save() methods by all ctools,
 * and it implements the execute() method.
 *
 * The class also provides a number of helper methods for ctools.
 ***************************************************************************/
class ctool : public GApplication  {
public:
    // Constructors and destructors
    ctool(const std::string& name, const std::string& version);
    ctool(const std::string& name, const std::string& version,
          int argc, char* argv[]);
    ctool(const ctool& app);
    virtual ~ctool(void);

    // Operators
    ctool& operator=(const ctool& app);

    // Pure virtual methods
    virtual void clear(void) = 0;
    virtual void run(void) = 0;
    virtual void save(void) = 0;

    // Public methods
    virtual void execute(void);

#ifdef SWIG
public:
#else
protected:
#endif
    // Protected methods
    void                  init_members(void);
    void                  copy_members(const ctool& app);
    void                  free_members(void);
    const bool&           read_ahead(void) const;
    const GTimeReference& time_reference(void) const;

    // Protected high-level setup methods
    void setup_observations(GObservations& obs, const bool& response = true,
                                                const bool& list = true,
                                                const bool& cube = true);
    void setup_models(GObservations& obs, const std::string& name = "");

    // Protected methods that create objects from user parameters
    GEbounds        create_ebounds(void);
    GSkyMap         create_map(const GObservations& obs);
    GCTAEventCube   create_cube(const GObservations& obs);
    GCTAObservation create_cta_obs(void);

    // Protected methods that check user parameters
    void            require_inobs(const std::string& method);
    void            require_inobs_nolist(const std::string& method);
    void            require_inobs_nocube(const std::string& method);

    // Protected methods for logging
    void            log_observations(const GChatter&      chatter,
                                     const GObservations& obs,
                                     const std::string&   what = "Observation");
    void            log_models(const GChatter&    chatter,
                               const GModels&     models,
                               const std::string& what = "Model");

    // Protected methods that extract user parameters
    GCTARoi           get_roi(void);

    // Protected support methods
    void              set_response(GObservations& obs);
    std::vector<bool> set_edisp(GObservations& obs, const bool& edisp) const;
    void              restore_edisp(GObservations& obs,
                                    const std::vector<bool>& edisp) const;
    void              set_obs_response(GCTAObservation* obs);
    void              set_obs_bounds(GObservations& obs);
    GObservations     get_observations(const bool& get_response = true);
    GSkyDir           get_mean_pointing(const GObservations& obs);
    size_t            get_current_rss(void);
    std::string       get_obs_header(const GObservation* obs) const;
    GEnergies         insert_energy_boundaries(const GEnergies&       energies,
                                               const GCTAObservation& obs);
    std::vector<bool> cube_layer_usage(const GEbounds& cube_ebounds,
                                       const GEbounds& list_ebounds) const;
    bool              is_valid_filename(const GFilename& filename) const;

    // Protected warning strings
    std::string     warn_too_few_energies(const GEnergies& energies) const;
    std::string     warn_xml_suffix(const GFilename& filename) const;

    // Protected members
    bool            m_read_ahead; //!< Read ahead output parameters

protected:
    // Protected methods
    void            provide_help(void) const;
    
    // Protected members
    bool            m_use_xml;  //!< Use XML file instead of FITS file for observations
    GTimeReference  m_cta_ref;  //!< CTA time reference
};


/***********************************************************************//**
 * @brief Signal whether parameters should be read ahead
 *
 * @return True if parameters should be read ahead
 *
 * Ahead reading of parameter is useful in case that some parameters are only
 * required in a late stage of the processing. By reading them ahead it is
 * assured that these parameters will be queried at the beginning of the
 * tool execution.
 ***************************************************************************/
inline
const bool& ctool::read_ahead(void) const
{
    return (m_read_ahead);
}


/***********************************************************************//**
 * @brief Return time reference
 *
 * @return Reference to time reference
 ***************************************************************************/
inline
const GTimeReference& ctool::time_reference(void) const
{
    return (m_cta_ref);
}


/***********************************************************************//**
 * @brief Check of a filename is valid
 *
 * @return True if filename is valid
 *
 * Checks whether a filename is not empty and not "NONE".
 ***************************************************************************/
inline
bool ctool::is_valid_filename(const GFilename& filename) const
{
    return (!filename.is_empty() &&
            (gammalib::tolower(filename.url()) != "none"));
}

#endif /* CTOOL_HPP */
