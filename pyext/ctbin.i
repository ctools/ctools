/***************************************************************************
 *                        ctbin - Event binning tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @file ctbin.i
 * @brief Event binning tool definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctbin.hpp"
%}


/***********************************************************************//**
 * @class ctbin
 *
 * @brief Event binning tool
 ***************************************************************************/
class ctbin : public ctobservation {
public:
    // Constructors and destructors
    ctbin(void);
    explicit ctbin(const GObservations& obs);
    ctbin(int argc, char *argv[]);
    ctbin(const ctbin& app);
    virtual ~ctbin(void);

    // Methods
    void                 clear(void);
    void                 run(void);
    int                  cubes(void) const;
    const GCTAEventCube& cube(const int& index = 0) const;
    void                 save(void);
    void                 publish(const std::string& name = "");
};


/***********************************************************************//**
 * @brief Event binning tool Python extension
 ***************************************************************************/
%extend ctbin {
    ctbin copy() {
        return (*self);
    }
}
