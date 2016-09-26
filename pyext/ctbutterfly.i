/***************************************************************************
 *                  ctbutterfly - butterfly calculation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Michael Mayer                                    *
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
 * @file ctbutterfly.i
 * @brief Butterfly calculation tool
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctbutterfly.hpp"
%}


/***********************************************************************//**
 * @class ctbutterfly
 *
 * @brief Butterfly calculation tool
 ***************************************************************************/
class ctbutterfly : public ctlikelihood {
public:
    // Constructors and destructors
    ctbutterfly(void);
    explicit ctbutterfly(const GObservations& obs);
    ctbutterfly(int argc, char *argv[]);
    ctbutterfly(const ctbutterfly& app);
    virtual ~ctbutterfly(void);

    // Methods
    void clear(void);
    void run(void);
    void save(void);
};


/***********************************************************************//**
 * @brief butterfly calculation tool Python extensions
 ***************************************************************************/
%extend ctbutterfly {
    ctbutterfly copy() {
        return (*self);
    }
}
