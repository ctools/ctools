/***************************************************************************
 *          ctedispcube - Energy dispersion cube generation tool           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Maria Haupt                                      *
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
 * @file ctedispcube.i
 * @brief Energy dispersion cube generation tool definition
 * @author Chia-Chun Lu
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctedispcube.hpp"
%}

/***********************************************************************//**
 * @class ctedispcube
 *
 * @brief Energy dispersion cube generation tool
 ***************************************************************************/
class ctedispcube : public ctobservation {

public:
    // Constructors and destructors
    ctedispcube(void);
    explicit ctedispcube(const GObservations& obs);
    ctedispcube(int argc, char *argv[]);
    ctedispcube(const ctedispcube& app);
    virtual ~ctedispcube(void);

    // Methods
    void                 clear(void);
    void                 run(void);
    void                 save(void);
    const GCTACubeEdisp& edispcube(void) const;
};

/***********************************************************************//**
 * @brief Energy dispersion cube generation tool Python extensions
 ***************************************************************************/
%extend ctedispcube {
    ctedispcube copy() {
        return (*self);
    }
}
