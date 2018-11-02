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
 * @file ctool_base.i
 * @brief [WHAT] tool definition
 * @author [AUTHOR]
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctool_base.hpp"
%}


/***********************************************************************//**
 * @class ctool_base
 *
 * @brief [WHAT] tool
 ***************************************************************************/
class ctool_base : public ctool {
public:
    // Constructors and destructors
    ctool_base(void);
    ctool_base(int argc, char *argv[]);
    ctool_base(const ctool_base& app);
    virtual ~ctool_base(void);

    // Methods
    void clear(void);
    void run(void);
    void save(void);
};


/***********************************************************************//**
 * @brief [WHAT] tool Python extension
 ***************************************************************************/
%extend ctool_base {
    ctool_base copy() {
        return (*self);
    }
}
