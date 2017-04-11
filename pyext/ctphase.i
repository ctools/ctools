/***************************************************************************
 *          ctphase - Append phase information to CTA events file          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Joshua Cardenzana                                *
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
 * @file ctphase.i
 * @brief Append phase information to CTA events file
 * @author Joshua Cardenzana
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctphase.hpp"
%}


/***********************************************************************//**
 * @class ctphase
 *
 * @brief Data selection tool
 ***************************************************************************/
class ctphase : public ctobservation {
public:
    // Constructors and destructors
    ctphase(void);
    explicit ctphase(const GObservations& obs);
    ctphase(int argc, char *argv[]);
    ctphase(const ctphase& app);
    virtual ~ctphase(void);
    
    // Methods
    void clear(void);
    void run(void);
    void save(void);
    void publish(const std::string& name = "");
};


/***********************************************************************//**
 * @brief Data selection tool Python extension
 ***************************************************************************/
%extend ctphase {
    ctphase copy() {
        return (*self);
    }
}
