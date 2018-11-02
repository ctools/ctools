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
 * @file ctool_like.i
 * @brief [WHAT] tool definition
 * @author [AUTHOR]
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctool_like.hpp"
%}


/***********************************************************************//**
 * @class ctool_like
 *
 * @brief [WHAT] tool
 ***************************************************************************/
class ctool_like : public ctlikelihood {
public:
    // Constructors and destructors
    ctool_like(void);
    explicit ctool_like(const GObservations& obs);
    ctool_like(int argc, char *argv[]);
    ctool_like(const ctool_like& app);
    virtual ~ctool_like(void);

    // Methods
    void clear(void);
    void run(void);
    void save(void);
};


/***********************************************************************//**
 * @brief [WHAT] tool Python extension
 ***************************************************************************/
%extend ctool_like {
    ctool_like copy() {
        return (*self);
    }
}
