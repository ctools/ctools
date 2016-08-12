/***************************************************************************
 *                 cterror - Parameter error calculation tool              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2016 by Florent Forest                              *
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
 * @file cterror.i
 * @brief Parameter error calculation tool interface definition
 * @author Florent Forest
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "cterror.hpp"
%}


/***********************************************************************//**
 * @class cterror
 *
 * @brief Parameter error calculation tool
 ***************************************************************************/
class cterror : public ctlikelihood {

public:
    // Constructors and destructors
    cterror(void);
    explicit cterror(const GObservations& obs);
    cterror(int argc, char *argv[]);
    cterror(const cterror& app);
    virtual ~cterror(void);

    // Methods
    void clear(void);
    void run(void);
    void save(void);
};


/***********************************************************************//**
 * @brief Parameter error calculation tool Python extensions
 ***************************************************************************/
%extend cterror {
    cterror copy() {
        return (*self);
    }
}
