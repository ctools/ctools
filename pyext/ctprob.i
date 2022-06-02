/***************************************************************************
 *           ctprob - Computes probability for a given model               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2022 by Leonardo Di Venere                          *
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
 * @file ctprob.i
 * @brief Computes probability for a given model
 * @author Leonardo Di Venere
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctprob.hpp"
%}


/***********************************************************************//**
 * @class ctprob
 *
 * @brief Data selection tool
 ***************************************************************************/
class ctprob : public ctobservation {
public:
    // Constructors and destructors
    ctprob(void);
    explicit ctprob(const GObservations& obs);
    ctprob(int argc, char *argv[]);
    ctprob(const ctprob& app);
    virtual ~ctprob(void);

    // Methods
    void clear(void);
    void process(void);
    void save(void);
    void publish(const std::string& name = "");
};


/***********************************************************************//**
 * @brief Data selection tool Python extension
 ***************************************************************************/
%extend ctprob {
    ctprob copy() {
        return (*self);
    }
}
