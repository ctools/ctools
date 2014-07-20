/***************************************************************************
 *                   ctexpcube - CTA exposure cube tool                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file ctexpcube.hpp
 * @brief CTA exposure cube tool definition
 * @author Juergen Knoedlseder
 */

#ifndef CTEXPCUBE_HPP
#define CTEXPCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTEXPCUBE_NAME    "ctexpcube"
#define CTEXPCUBE_VERSION "00-01-00"


/***********************************************************************//**
 * @class ctexpcube
 *
 * @brief CTA exposure cube tool
 ***************************************************************************/
class ctexpcube : public GApplication  {

public:
    // Constructors and destructors
    ctexpcube(void);
    explicit ctexpcube(const GObservations& obs);
    ctexpcube(int argc, char *argv[]);
    ctexpcube(const ctexpcube& app);
    virtual ~ctexpcube(void);

    // Operators
    ctexpcube& operator=(const ctexpcube& app);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           get_obs(void);
    void           set_response(void);
    void           get_ebounds(void);
    void           set_from_cntmap(const std::string& filename);

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const ctexpcube& app);
    void           free_members(void);

    // User parameters
    std::string   m_outfile;     //!< Output exposure cube file
    bool          m_apply_edisp; //!< Apply energy dispersion?

    // Protected members
    bool          m_read_ahead;  //!< Read ahead parameters
    GObservations m_obs;         //!< Observation container
    GCTAExposure  m_expcube;     //!< Exposure cube
    GEbounds      m_ebounds;     //!< Energy boundaries
};

#endif /* CTEXPCUBE_HPP */
