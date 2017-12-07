/***************************************************************************
 *                   ctulimit - Upper limit calculation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2016 by Michael Mayer                               *
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
 * @file ctulimit
 * @brief Upper limit calculation tool
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctulimit.hpp"
%}


/***********************************************************************//**
 * @class ctulimit
 *
 * @brief Upper limit calculation tool
 ***************************************************************************/
class ctulimit : public ctlikelihood {
public:
    // Constructors and destructors
    ctulimit(void);
    explicit ctulimit(const GObservations& obs);
    ctulimit(int argc, char *argv[]);
    ctulimit(const ctulimit& app);
    virtual ~ctulimit(void);
    
    // Methods
    void          clear(void);
    void          run(void);
    void          save(void);
    const double& diff_ulimit(void) const;
    const double& flux_ulimit(void) const;
    const double& eflux_ulimit(void) const;
};


/***********************************************************************//**
 * @brief Upper limit calculation tool Python extensions
 ***************************************************************************/
%extend ctulimit {
    ctulimit copy() {
        return (*self);
    }
}
