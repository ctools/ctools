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
#define SWIG
#include "ctobservation.hpp"


/***********************************************************************//**
 * @class csobservation
 *
 * @brief csobservation base class
 *
 * This class is a non-abstract C++ implementation of the ctobservation base
 * class. It serves as Python base class since Python does not know about
 * abstract classes.
 ***************************************************************************/
class csobservation : public ctobservation  {
public:
    // Constructors and destructors
    csobservation(const std::string& name, const std::string& version) :
                  ctobservation(name, version) {}
    csobservation(const std::string& name, const std::string& version,
                  const GObservations& obs) :
                  ctobservation(name, version, obs) {}
    csobservation(const std::string& name, const std::string& version,
                  int argc, char* argv[]) :
                  ctobservation(name, version, argc, argv){}
    csobservation(const csobservation& app) : ctobservation(app) {}
    virtual ~csobservation(void) {}

    // Dummy methods (implementation makes this class non-abstract)
    virtual void clear(void) {}
    virtual void run(void) {}
    virtual void save(void) {}
};
%}

// Include (int ARGC, char **ARGV) typemap to allow passing command line
// arguments to GApplication constructor
%include "argcargv.i"


/***********************************************************************//**
 * @class ctobservation
 *
 * @brief Base class for observation tools
 ***************************************************************************/
class ctobservation : public ctool {

public:
    // Constructors and destructors
    ctobservation(const std::string& name, const std::string& version);
    ctobservation(const std::string& name, const std::string& version,
                  const GObservations& obs);
    ctobservation(const std::string& name, const std::string& version,
                  int ARGC, char **ARGV);
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
 * @class csobservation
 *
 * @brief Base class for observation scripts
 *
 * This is the base class from which observation scripts should derive.
 ***************************************************************************/
class csobservation : public ctobservation  {
public:        
    // Constructors and destructors
    csobservation(const std::string& name, const std::string& version);
    csobservation(const std::string& name, const std::string& version,
                  const GObservations& obs);
    csobservation(const std::string& name, const std::string& version,
                  int ARGC, char **ARGV);
    csobservation(const csobservation& app);
    virtual ~csobservation(void);

    // Methods
    virtual void clear(void);
    virtual void run(void);
    virtual void save(void);
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

# Initialise application by calling the appropriate base class constructor.
# The function supports either an observation container, or an argument
# list or no argument as "argv" parameter. The function also writes the
# header in the log file and switches the date on for logging.
def _init_csobservation(self, argv):
    if len(argv) > 0 and isinstance(argv[0],gammalib.GObservations):
        csobservation.__init__(self, self._name, self._version, argv[0])
    elif len(argv) > 0:
        csobservation.__init__(self, self._name, self._version, *argv)
    else:
        csobservation.__init__(self, self._name, self._version)
    # Set logger properties
    self._log_header()
    self._log.date(True)
csobservation._init_csobservation = _init_csobservation

# Define an iterator over all observations
def _unbinned_observations(self):
    obs = self._first_unbinned_observation()
    while obs != None:
        yield obs
        obs = self._next_unbinned_observation()
ctool._unbinned_observations   = _unbinned_observations
cscript._unbinned_observations = _unbinned_observations
%}
