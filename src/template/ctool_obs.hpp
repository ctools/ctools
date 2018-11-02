/***************************************************************************
 *                        ctool_obs - [WHAT] tool                          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file ctool_obs.hpp
 * @brief [WHAT] tool definition
 * @author [AUTHOR]
 */

#ifndef CTOOL_OBS_HPP
#define CTOOL_OBS_HPP

/* __ Includes ___________________________________________________________ */
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTOOL_OBS_NAME "ctool_obs"


/***********************************************************************//**
 * @class ctool_obs
 *
 * @brief [WHAT] tool
 *
 * @todo Add tool description.
 ***************************************************************************/
class ctool_obs : public ctobservation {
public:
    // Constructors and destructors
    ctool_obs(void);
    explicit ctool_obs(const GObservations& obs);
    ctool_obs(int argc, char *argv[]);
    ctool_obs(const ctool_obs& app);
    virtual ~ctool_obs(void);

    // Operators
    ctool_obs& operator=(const ctool_obs& app);

    // Methods
    void clear(void);
    void run(void);
    void save(void);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctool_obs& app);
    void free_members(void);
    void get_parameters(void);

    // Protected members
    // TODO: Add any data members that are necessary
};

#endif /* CTOOL_OBS_HPP */
