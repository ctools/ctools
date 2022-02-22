/***************************************************************************
 *   ctfindvar - search time variability tool                              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018-2022 by Simon Bonnefoy                              *
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
 * @file ctfindvar.i
 * @brief search time variability tool definition
 * @author Simon Bonnefoy
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctfindvar.hpp"
%}


/***********************************************************************//**
 * @class ctfindvar
 *
 * @brief search time variability tool
 ***************************************************************************/
class ctfindvar : public ctobservation {
public:
    // Constructors and destructors
    ctfindvar(void);
    explicit ctfindvar(const GObservations& obs);
    ctfindvar(int argc, char *argv[]);
    ctfindvar(const ctfindvar& app);
    virtual ~ctfindvar(void);

    // Methods
    void clear(void);
    void process(void);
    void save(void);
};


/***********************************************************************//**
 * @brief search time variability tool Python extension
 ***************************************************************************/
%extend ctfindvar {
    ctfindvar copy() {
        return (*self);
    }
}
