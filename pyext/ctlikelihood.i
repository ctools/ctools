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
#define SWIG
#include "ctlikelihood.hpp"


/***********************************************************************//**
 * @class cslikelihood
 *
 * @brief cslikelihood base class
 *
 * This class is a non-abstract C++ implementation of the ctlikelihood base
 * class. It serves as Python base class since Python does not know about
 * abstract classes.
 ***************************************************************************/
class cslikelihood : public ctlikelihood  {
public:
    // Constructors and destructors
    cslikelihood(const std::string& name, const std::string& version) :
                 ctlikelihood(name, version) {}
    cslikelihood(const std::string& name, const std::string& version,
                 const GObservations& obs) :
                 ctlikelihood(name, version, obs) {}
    cslikelihood(const std::string& name, const std::string& version,
                 int argc, char* argv[]) :
                 ctlikelihood(name, version, argc, argv){}
    cslikelihood(const cslikelihood& app) : ctlikelihood(app) {}
    virtual ~cslikelihood(void) {}

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

    // Make methods private in Python by prepending an underscore
    %rename(_evaluate) evaluate;

    // Protected methods
    double evaluate(GModelPar& par, const double& value);
};


/***********************************************************************//**
 * @class cslikelihood
 *
 * @brief Base class for likelihood scripts
 *
 * This is the base class from which likelihood scripts should derive.
 ***************************************************************************/
class cslikelihood : public ctlikelihood  {
public:        
    // Constructors and destructors
    cslikelihood(const std::string& name, const std::string& version);
    cslikelihood(const std::string& name, const std::string& version,
                 const GObservations& obs);
    cslikelihood(const std::string& name, const std::string& version,
                 int ARGC, char **ARGV);
    cslikelihood(const cslikelihood& app);
    virtual ~cslikelihood(void);

    // Methods
    virtual void clear(void);
    virtual void run(void);
    virtual void save(void);
};


/***********************************************************************//**
 * @brief Likelihood tool base class C++ extensions
 ***************************************************************************/
%extend ctlikelihood {
}


/***********************************************************************//**
 * @brief Likelihood tool base class Python extensions
 ***************************************************************************/
%pythoncode %{

# Initialise application by calling the appropriate base class constructor.
# The function supports either an observation container, or an argument
# list or no argument as "argv" parameter. The function also writes the
# header in the log file and switches the date on for logging.
def _init_cslikelihood(self, argv):
    if len(argv) > 0 and isinstance(argv[0],gammalib.GObservations):
        cslikelihood.__init__(self, self._name, self._version, argv[0])
    elif len(argv) > 0:
        cslikelihood.__init__(self, self._name, self._version, *argv)
    else:
        cslikelihood.__init__(self, self._name, self._version)
    # Set logger properties
    self._log_header()
    self._log.date(True)
cslikelihood._init_cslikelihood = _init_cslikelihood
%}
