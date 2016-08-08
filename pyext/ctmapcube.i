/***************************************************************************
 *                  ctmapcube - Map cube generation tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file ctmapcube.i
 * @brief Map cube generation tool definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctmapcube.hpp"
%}


/***********************************************************************//**
 * @class ctmapcube
 *
 * @brief Map cube generation tool
 ***************************************************************************/
class ctmapcube : public ctool {

public:
    // Constructors and destructors
    ctmapcube(void);
    ctmapcube(int argc, char *argv[]);
    ctmapcube(const ctmapcube& app);
    virtual ~ctmapcube(void);

    // Methods
    void                            clear(void);
    void                            run(void);
    void                            save(void);
    void                            publish(const std::string& name = "");
    const GModelSpatialDiffuseCube& mapcube(void) const;
    void                            models(const GModels& models);
};


/***********************************************************************//**
 * @brief Map cube generation tool Python extension
 ***************************************************************************/
%extend ctmapcube {
    ctmapcube copy() {
        return (*self);
    }
}
