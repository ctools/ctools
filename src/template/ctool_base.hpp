/***************************************************************************
 *                        ctool_base - [WHAT] tool                         *
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
 * @file ctool_base.hpp
 * @brief [WHAT] tool definition
 * @author [AUTHOR]
 */

#ifndef CTOOL_BASE_HPP
#define CTOOL_BASE_HPP

/* __ Includes ___________________________________________________________ */
#include "ctool.hpp"

/* __Definitions _________________________________________________________ */
#define CTOOL_BASE_NAME "ctool_base"


/***********************************************************************//**
 * @class ctool_base
 *
 * @brief [WHAT] tool
 *
 * @todo Add tool description.
 ***************************************************************************/
class ctool_base : public ctool {
public:
    // Constructors and destructors
    ctool_base(void);
    ctool_base(int argc, char *argv[]);
    ctool_base(const ctool_base& app);
    virtual ~ctool_base(void);

    // Operators
    ctool_base& operator=(const ctool_base& app);

    // Methods
    void clear(void);
    void process(void);
    void save(void);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctool_base& app);
    void free_members(void);
    void get_parameters(void);

    // Protected members
    // TODO: Add any data members that are necessary
};

#endif /* CTOOL_BASE_HPP */
