/***************************************************************************
 *                  ctmodel - Model cube generation tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * @file ctmodel.i
 * @brief Model cube generation tool definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctmodel.hpp"
%}


/***********************************************************************//**
 * @class ctmodel
 *
 * @brief Model cube generation tool
 ***************************************************************************/
class ctmodel : public ctobservation {
public:
    // Constructors and destructors
    ctmodel(void);
    explicit ctmodel(const GObservations& obs);
    ctmodel(int argc, char *argv[]);
    ctmodel(const ctmodel& app);
    virtual ~ctmodel(void);

    // Methods
    void                 clear(void);
    void                 run(void);
    void                 save(void);
    void                 publish(const std::string& name = "");
    const GCTAEventCube& cube(void) const;
    void                 cube(const GCTAEventCube& cube);
    void                 models(const GModels& models);
};


/***********************************************************************//**
 * @brief Model cube generation tool Python extension
 ***************************************************************************/
%extend ctmodel {
    ctmodel copy() {
        return (*self);
    }
}
