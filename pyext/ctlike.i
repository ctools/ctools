/***************************************************************************
 *                ctlike - Maximum likelihood fitting tool                 *
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
 * @file ctlike.i
 * @brief Maximum likelihood fitting tool definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctlike.hpp"
%}


/***********************************************************************//**
 * @class ctlike
 *
 * @brief Maximum likelihood fitting tool
 ***************************************************************************/
class ctlike : public ctlikelihood {
public:
    // Constructors and destructors
    ctlike(void);
    explicit ctlike(const GObservations& obs);
    ctlike(int argc, char *argv[]);
    ctlike(const ctlike& app);
    virtual ~ctlike(void);

    // Methods
    void clear(void);
    void run(void);
    void save(void);
};


/***********************************************************************//**
 * @brief Maximum likelihood fitting tool Python extension
 ***************************************************************************/
%extend ctlike {
    ctlike copy() {
        return (*self);
    }
}
