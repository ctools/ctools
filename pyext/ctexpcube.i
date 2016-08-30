/***************************************************************************
 *                 ctexpcube - Exposure cube generation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Juergen Knoedlseder                         *
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
 * @file ctexpcube.i
 * @brief Exposure cube generation tool definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctexpcube.hpp"
%}


/***********************************************************************//**
 * @class ctexpcube
 *
 * @brief Exposure cube generation tool
 ***************************************************************************/
class ctexpcube : public ctobservation {
public:
    // Constructors and destructors
    ctexpcube(void);
    explicit ctexpcube(const GObservations& obs);
    ctexpcube(int argc, char *argv[]);
    ctexpcube(const ctexpcube& app);
    virtual ~ctexpcube(void);

    // Methods
    void                    clear(void);
    void                    run(void);
    void                    save(void);
    void                    publish(const std::string& name = "");
    const GCTACubeExposure& expcube(void) const;
};


/***********************************************************************//**
 * @brief Exposure cube generation tool Python extension
 ***************************************************************************/
%extend ctexpcube {
    ctexpcube copy() {
        return (*self);
    }
}
