/***************************************************************************
 *                           xxx - [WHAT] tool                             *
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
 * @file xxx.hpp
 * @brief [WHAT] tool definition
 * @author [AUTHOR]
 */

#ifndef XXX_HPP
#define XXX_HPP

/* __ Includes ___________________________________________________________ */
#include "ctool.hpp"

/* __Definitions _________________________________________________________ */
#define XXX_NAME "xxx"


/***********************************************************************//**
 * @class xxx
 *
 * @brief [WHAT] tool
 *
 * @todo Add tool description.
 ***************************************************************************/
class xxx : public ctool {
public:
    // Constructors and destructors
    xxx(void);
    xxx(int argc, char *argv[]);
    xxx(const xxx& app);
    virtual ~xxx(void);

    // Operators
    xxx& operator=(const xxx& app);

    // Methods
    void clear(void);
    void run(void);
    void save(void);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const xxx& app);
    void free_members(void);
    void get_parameters(void);

    // Protected members
    // TODO: Add any data members that are necessary
};

#endif /* XXX_HPP */
