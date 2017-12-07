/***************************************************************************
 *                    cttsmap - TS map calculation tool                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Michael Mayer                               *
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
 * @file cttsmap.i
 * @brief TS map calculation tool interface definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "cttsmap.hpp"
%}


/***********************************************************************//**
 * @class cttsmap
 *
 * @brief TS map calculation tool
 ***************************************************************************/
class cttsmap : public ctlikelihood {

public:
    // Constructors and destructors
    cttsmap(void);
    explicit cttsmap(const GObservations& obs);
    cttsmap(int argc, char *argv[]);
    cttsmap(const cttsmap& app);
    virtual ~cttsmap(void);
    
    // Methods
    void           clear(void);
    void           run(void);
    void           save(void);
    void           publish(const std::string& name = "");
    const GSkyMap& tsmap(void) const;
};


/***********************************************************************//**
 * @brief TS map calculation tool Python extensions
 ***************************************************************************/
%extend cttsmap {
    cttsmap copy() {
        return (*self);
    }
}
