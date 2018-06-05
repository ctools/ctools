/***************************************************************************
 *                      ctselect - Data selection tool                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Jurgen Knodlseder                           *
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
 * @file ctselect.i
 * @brief Data selection tool definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctselect.hpp"
%}


/***********************************************************************//**
 * @class ctselect
 *
 * @brief Data selection tool
 ***************************************************************************/
class ctselect : public ctobservation {
public:
    // Constructors and destructors
    ctselect(void);
    explicit ctselect(const GObservations& obs);
    ctselect(int argc, char *argv[]);
    ctselect(const ctselect& app);
    virtual ~ctselect(void);

    // Methods
    void clear(void);
    void run(void);
    void save(void);
    void publish(const std::string& name = "");
};


/***********************************************************************//**
 * @brief Data selection tool Python extension
 ***************************************************************************/
%extend ctselect {
    ctselect copy() {
        return (*self);
    }
}
