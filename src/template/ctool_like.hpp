/***************************************************************************
 *                       ctool_like - [WHAT] tool                          *
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
 * @file ctool_like.hpp
 * @brief [WHAT] tool definition
 * @author [AUTHOR]
 */

#ifndef CTOOL_LIKE_HPP
#define CTOOL_LIKE_HPP

/* __ Includes ___________________________________________________________ */
#include "ctlikelihood.hpp"

/* __Definitions _________________________________________________________ */
#define CTOOL_LIKE_NAME "ctool_like"


/***********************************************************************//**
 * @class ctool_like
 *
 * @brief [WHAT] tool
 *
 * @todo Add tool description.
 ***************************************************************************/
class ctool_like : public ctlikelihood {
public:
    // Constructors and destructors
    ctool_like(void);
    explicit ctool_like(const GObservations& obs);
    ctool_like(int argc, char *argv[]);
    ctool_like(const ctool_like& app);
    virtual ~ctool_like(void);

    // Operators
    ctool_like& operator=(const ctool_like& app);

    // Methods
    void clear(void);
    void process(void);
    void save(void);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctool_like& app);
    void free_members(void);
    void get_parameters(void);

    // Protected members
    // TODO: Add any data members that are necessary
};

#endif /* CTOOL_LIKE_HPP */
