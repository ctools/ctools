/***************************************************************************
 *                  ctobssim - Observation simulator tool                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2016 by Juergen Knoedlseder                         *
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
 * @file ctobssim.i
 * @brief Observation simulator tool definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctobssim.hpp"
%}


/***********************************************************************//**
 * @class ctobssim
 *
 * @brief Observation simulator tool
 ***************************************************************************/
class ctobssim : public ctobservation {
public:
    // Constructors and destructors
    ctobssim(void);
    explicit ctobssim(const GObservations& obs);
    ctobssim(int argc, char *argv[]);
    ctobssim(const ctobssim& app);
    virtual ~ctobssim(void);

    // Methods
    void          clear(void);
    void          run(void);
    void          save(void);
    void          publish(const std::string& name = "");
    const double& max_rate(void) const;
    void          max_rate(const double& max_rate);
};


/***********************************************************************//**
 * @brief Observation simulation tool Python extension
 ***************************************************************************/
%extend ctobssim {
    ctobssim copy() {
        return (*self);
    }
}
