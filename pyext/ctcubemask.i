/***************************************************************************
 *                      ctcubemask - Cube filter tool                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 Chia-Chun Lu                                   *
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
 * @file ctcubemask.i
 * @brief Cube filter tool definition
 * @author Chia-Chun Lu
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctcubemask.hpp"
%}


/***********************************************************************//**
 * @class ctcubemask
 *
 * @brief Cube filter tool
 ***************************************************************************/
class ctcubemask : public ctobservation {
public:
    // Constructors and destructors
    ctcubemask(void);
    explicit ctcubemask(const GObservations& obs);
    ctcubemask(int argc, char *argv[]);
    ctcubemask(const ctcubemask& app);
    virtual ~ctcubemask(void);

    // Methods
    void clear(void);
    void run(void);
    void save(void);
    void publish(const std::string& name = "");
};


/***********************************************************************//**
 * @brief Cube filter tool Python extension
 ***************************************************************************/
%extend ctcubemask {
    ctcubemask copy() {
        return (*self);
    }
}
