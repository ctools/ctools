/***************************************************************************
 *              ctlikelihood - Base class for likelihood tools             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file ctlikelihood.i
 * @brief Likelihood tool base class interface definition
 * @author Juergen Knoedlseder
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctlikelihood.hpp"
%}


/***********************************************************************//**
 * @class ctlikelihood
 *
 * @brief Base class for likelihood tools
 ***************************************************************************/
class ctlikelihood : public ctobservation {

public:
    // Constructors and destructors
    ctlikelihood(const std::string& name, const std::string& version);
    ctlikelihood(const std::string& name, const std::string& version,
                 const GObservations& obs);
    ctlikelihood(const std::string& name, const std::string& version,
                 int argc, char* argv[]);
    ctlikelihood(const ctlikelihood& app);
    virtual ~ctlikelihood(void);

    // Pure virtual methods
    virtual void clear(void) = 0;
    virtual void run(void) = 0;
    virtual void save(void) = 0;

    // Methods
    const GOptimizer* opt(void) const;
};


/***********************************************************************//**
 * @brief Likelihood tool Python extensions
 ***************************************************************************/
%extend ctlikelihood {
}
