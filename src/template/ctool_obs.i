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
 * @file xxx.i
 * @brief [WHAT] tool definition
 * @author [AUTHOR]
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "xxx.hpp"
%}


/***********************************************************************//**
 * @class xxx
 *
 * @brief [WHAT] tool
 ***************************************************************************/
class xxx : public ctobservation {
public:
    // Constructors and destructors
    xxx(void);
    explicit xxx(const GObservations& obs);
    xxx(int argc, char *argv[]);
    xxx(const xxx& app);
    virtual ~xxx(void);

    // Methods
    void clear(void);
    void run(void);
    void save(void);
};


/***********************************************************************//**
 * @brief [WHAT] tool Python extension
 ***************************************************************************/
%extend xxx {
    xxx copy() {
        return (*self);
    }
}
