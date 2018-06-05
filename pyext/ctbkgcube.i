/***************************************************************************
 *               ctbkgcube - Background cube generation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Chia-Chun Lu                                *
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
 * @file ctbkgcube.i
 * @brief Background cube generation tool definition
 * @author Chia-Chun Lu
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctbkgcube.hpp"
%}

/***********************************************************************//**
 * @class ctbkgcube
 *
 * @brief Background cube generation tool
 ***************************************************************************/
class ctbkgcube : public ctobservation  {
public:
    // Constructors and destructors
    ctbkgcube(void);
    explicit ctbkgcube(const GObservations& obs);
    ctbkgcube(int argc, char *argv[]);
    ctbkgcube(const ctbkgcube& app);
    virtual ~ctbkgcube(void);

    // Methods
    void                      clear(void);
    void                      run(void);
    void                      save(void);
    void                      publish(const std::string& name = "");
    const GCTACubeBackground& bkgcube(void) const;
    const GModels&            models(void) const;
};

/***********************************************************************//**
 * @brief Background cube generation tool Python extension
 ***************************************************************************/
%extend ctbkgcube {
    ctbkgcube copy() {
        return (*self);
    }
}
