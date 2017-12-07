/***************************************************************************
 *                  ctpsfcube - PSF cube generation tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Chia-Chun Lu                                *
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
 * @file ctpsfcube.i
 * @brief PSF cube generation tool definition
 * @author Chia-Chun Lu
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctpsfcube.hpp"
%}

/***********************************************************************//**
 * @class ctpsfcube
 *
 * @brief PSF cube generation tool
 ***************************************************************************/
class ctpsfcube : public ctobservation {

public:
    // Constructors and destructors
    ctpsfcube(void);
    explicit ctpsfcube(const GObservations& obs);
    ctpsfcube(int argc, char *argv[]);
    ctpsfcube(const ctpsfcube& app);
    virtual ~ctpsfcube(void);

    // Methods
    void               clear(void);
    void               run(void);
    void               save(void);
    const GCTACubePsf& psfcube(void) const;
};

/***********************************************************************//**
 * @brief PSF cube generation tool Python extension
 ***************************************************************************/
%extend ctpsfcube {
    ctpsfcube copy() {
        return (*self);
    }
}
