/***************************************************************************
 *                  ctmapcube - Map cube generation tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2017 by Juergen Knoedlseder                         *
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
 * @file ctmapcube.hpp
 * @brief Map cube generation tool definition
 * @author Juergen Knoedlseder
 */

#ifndef CTMAPCUBE_HPP
#define CTMAPCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctool.hpp"

/* __Definitions _________________________________________________________ */
#define CTMAPCUBE_NAME    "ctmapcube"
#define CTMAPCUBE_VERSION "1.3.0"


/***********************************************************************//**
 * @class ctmapcube
 *
 * @brief Map cube generation tool
 ***************************************************************************/
class ctmapcube : public ctool {

public:
    // Constructors and destructors
    ctmapcube(void);
    ctmapcube(int argc, char *argv[]);
    ctmapcube(const ctmapcube& app);
    virtual ~ctmapcube(void);

    // Operators
    ctmapcube& operator=(const ctmapcube& app);

    // Methods
    void                            clear(void);
    void                            run(void);
    void                            save(void);
    void                            publish(const std::string& name = "");
    const GModelSpatialDiffuseCube& mapcube(void) const;
    void                            models(const GModels& models);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctmapcube& app);
    void free_members(void);
    void get_parameters(void);
    void create_cube(void);
    void add_model(GModelSky* model);
    void add_ptsrc_model(GModelSky* model);
    void get_bounding_circle(GModelSky* model,
                             GSkyDir*   dir,
                             double*    radius) const;

    // User parameters
    GFilename m_outcube;                //!< Output map cube filename
    double    m_ptsrcsig;               //!< Point source sigma (arcmin)
    bool      m_publish;                //!< Publish map cube?
    GChatter  m_chatter;                //!< Chattiness

    // Protected members
    GModels                  m_models;  //!< Model container
    GModelSpatialDiffuseCube m_cube;    //!< Map cube
};


/***********************************************************************//**
 * @brief Return map cube
 *
 * @return Map cube
 ***************************************************************************/
inline
const GModelSpatialDiffuseCube& ctmapcube::mapcube(void) const
{
    return (m_cube);
}


/***********************************************************************//**
 * @brief Set models
 *
 * @param[in] models Model container.
 *
 * Set model container to be used for map cube generation.
 ***************************************************************************/
inline
void ctmapcube::models(const GModels& models)
{
    m_models = models;
    return;
}

#endif /* CTMAPCUBE_HPP */
