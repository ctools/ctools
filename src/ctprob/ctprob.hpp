/***************************************************************************
 *                  ctmodel - Model cube generation tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * @file ctmodel.hpp
 * @brief Model cube generation tool definition
 * @author Juergen Knoedlseder
 */

#ifndef CTMODEL_HPP
#define CTMODEL_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTMODEL_NAME    "ctmodel"
#define CTMODEL_VERSION "1.1.0"


/***********************************************************************//**
 * @class ctmodel
 *
 * @brief Model cube generation tool
 *
 * This class creates counts model map(s). The definition of the counts model
 * can be taken from a predefined observation container, from the counts maps
 * found in an observation definition XML file, from an individual counts
 * map, or from the parameter definition.
 *
 * Results are stored in an observation container that can be written to disk
 * in form of FITS files (model maps) and an updated observation definition
 * XML file.
 ***************************************************************************/
class ctmodel : public ctobservation {

public:
    // Constructors and destructors
    ctmodel(void);
    explicit ctmodel(const GObservations& obs);
    ctmodel(int argc, char *argv[]);
    ctmodel(const ctmodel& app);
    virtual ~ctmodel(void);

    // Operators
    ctmodel& operator=(const ctmodel& app);

    // Methods
    void                 clear(void);
    void                 run(void);
    void                 save(void);
    void                 publish(const std::string& name = "");
    const GCTAEventCube& cube(void) const;
    void                 cube(const GCTAEventCube& cube);
    void                 models(const GModels& models);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctmodel& app);
    void free_members(void);
    void get_parameters(void);
    void get_obs(void);
    void fill_cube(const GCTAObservation* obs);
    bool has_cube(void) const;

    // User parameters
    GFilename m_outcube;      //!< Output model cube
    bool      m_apply_edisp;  //!< Apply energy dispersion?
    bool      m_publish;      //!< Publish model cube?
    GChatter  m_chatter;      //!< Chattiness

    // Protected members
    GCTAEventCube m_cube;        //!< Model cube
    GGti          m_gti;         //!< Model cube GTIs
    bool          m_has_cube;    //!< Signal if cube has been set or loaded
    bool          m_append_cube; //!< Signal that cube should be appended
    bool          m_binned;      //!< Signals that we are in binned mode
};


/***********************************************************************//**
 * @brief Return model cube
 *
 * @return Reference to model cube.
 ***************************************************************************/
inline
const GCTAEventCube& ctmodel::cube(void) const
{
    return m_cube;
}


/***********************************************************************//**
 * @brief Signal if cube has been set or loaded
 *
 * @return True if cube has been set or loaded.
 ***************************************************************************/
inline
bool ctmodel::has_cube(void) const
{
    return m_has_cube;
}


/***********************************************************************//**
 * @brief Set models
 *
 * @param[in] models Model container.
 *
 * Set model container that should be used for model generation.
 ***************************************************************************/
inline
void ctmodel::models(const GModels& models)
{
    m_obs.models(models);
    return;
}

#endif /* CTMODEL_HPP */
