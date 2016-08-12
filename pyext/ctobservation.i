/***************************************************************************
 *             ctobservation - Base class for observation tools            *
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
 * @file ctobservation.i
 * @brief Observation tool base class interface definition
 * @author Juergen Knoedlseder
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctobservation.hpp"
%}


/***********************************************************************//**
 * @class ctobservation
 *
 * @brief Base class for likelihood tools
 ***************************************************************************/
class ctobservation : public ctool {

public:
    // Constructors and destructors
    ctobservation(const std::string& name, const std::string& version);
    ctobservation(const std::string& name, const std::string& version,
                  const GObservations& obs);
    ctobservation(const std::string& name, const std::string& version,
                  int argc, char* argv[]);
    ctobservation(const ctobservation& app);
    virtual ~ctobservation(void);

    // Pure virtual methods
    virtual void clear(void) = 0;
    virtual void run(void) = 0;
    virtual void save(void) = 0;

    // Public methods
    const GObservations& obs(void) const;

    // Make methods private in Python by prepending an underscore
    %rename(_first_unbinned_observation) first_unbinned_observation;
    %rename(_next_unbinned_observation)  next_unbinned_observation;

    // Protected methods
    GCTAObservation* first_unbinned_observation(void);
    GCTAObservation* next_unbinned_observation(void);
};


/***********************************************************************//**
 * @brief Observation tool base class C++ extensions
 ***************************************************************************/
%extend ctobservation {
}


/***********************************************************************//**
 * @brief Observation tool base class Python extensions
 ***************************************************************************/
%pythoncode %{
def _unbinned_observations(self):
    obs = self._first_unbinned_observation()
    while obs != None:
        yield obs
        obs = self._next_unbinned_observation()
ctool._unbinned_observations   = _unbinned_observations
cscript._unbinned_observations = _unbinned_observations
%}
