/***************************************************************************
 *             ctobservation - Base class for observation tools            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2022 by Juergen Knoedlseder                         *
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
    csobservation(const std::string& name,
                  const std::string& version) :
                  ctobservation(name, version) {}
    csobservation(const std::string&      name,
                  const std::string&      version,
                  const GApplicationPars& pars) :
                  ctobservation(name, version, pars) {}
    csobservation(const std::string& name,
                  const std::string& version,
                  int                argc,
                  char*              argv[]) :
                  ctobservation(name, version, argc, argv) {}
    csobservation(const std::string&   name,
                  const std::string&   version,
                  const GObservations& obs) :
                  ctobservation(name, version, obs) {}
    csobservation(const csobservation& app) : ctobservation(app) {}
    virtual ~csobservation(void) {}

    // Dummy methods (implementation makes this class non-abstract)
    virtual void clear(void) {}
    virtual void process(void) {}
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
    ctobservation(const std::string& name,
                  const std::string& version);
    ctobservation(const std::string&      name,
                  const std::string&      version,
                  const GApplicationPars& pars);
    ctobservation(const std::string& name,
                  const std::string& version,
                  int                ARGC,
                  char               **ARGV);
    ctobservation(const std::string&   name,
                  const std::string&   version,
                  const GObservations& obs);
    ctobservation(const ctobservation& app);
    virtual ~ctobservation(void);

    // Pure virtual methods
    virtual void clear(void) = 0;
    virtual void process(void) = 0;
    virtual void save(void) = 0;

    // Public methods
    void                 obs(const GObservations& obs);
    const GObservations& obs(void) const;

    // Make methods private in Python by prepending an underscore
    %rename(_first_unbinned_observation) first_unbinned_observation;
    %rename(_next_unbinned_observation)  next_unbinned_observation;
    %rename(_read_ogip_keywords)         read_ogip_keywords;
    %rename(_write_ogip_keywords)        write_ogip_keywords;
    %rename(_set_obs_statistic)          set_obs_statistic;
    %rename(_set_obs_bounds)             set_obs_bounds;
    %rename(_save_events_fits)           save_events_fits;
    %rename(_save_events_xml)            save_events_xml;

    // Protected methods
    GCTAObservation* first_unbinned_observation(void);
    GCTAObservation* next_unbinned_observation(void);
    void             read_ogip_keywords(GFitsHDU* hdu) const;
    void             write_ogip_keywords(GFitsHDU* hdu) const;
    void             set_obs_statistic(const std::string& statistic);
    void             set_obs_bounds();
    void             save_events_fits(void);
    void             save_events_xml(void);
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
    csobservation(const std::string& name,
                  const std::string& version);
    csobservation(const std::string&      name,
                  const std::string&      version,
                  const GApplicationPars& pars);
    csobservation(const std::string& name,
                  const std::string& version,
                  int                ARGC,
                  char               **ARGV);
    csobservation(const std::string&   name,
                  const std::string&   version,
                  const GObservations& obs);
    csobservation(const csobservation& app);
    virtual ~csobservation(void);

    // Methods
    virtual void clear(void);
    virtual void process(void);
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
%extend csobservation {
%pythoncode {
    def __getstate__(self):
        state = (ctool.__getstate__(self), self.obs())
        return state
    def __setstate__(self, state):
        ctool.__setstate__(self, state[0])
        self.obs(state[1])
}
}


/***********************************************************************//**
 * @brief Observation tool base class Python extensions
 ***************************************************************************/
%pythoncode %{

# Initialise application by calling the appropriate base class constructor.
# The function supports either an observation container, or an argument
# list or no argument as "argv" parameter. The function also writes the
# header in the log file and switches the date on for logging.
def _init_csobservation(self, name, version, argv):
    if len(argv) > 0 and isinstance(argv[0],gammalib.GObservations):
        csobservation.__init__(self, name, version, argv[0])
    elif len(argv) > 0:
        if len(argv) == 3 and argv[0] == name and argv[1] == version:
            csobservation.__init__(self, name, version, argv[2])
        else:
            csobservation.__init__(self, name, version, *argv)
    else:
        csobservation.__init__(self, name, version)
    # Set logger properties
    self._log.date(True)
csobservation._init_csobservation = _init_csobservation

# Define an iterator over all unbinned observations
def _unbinned_observations(self):
    obs = self._first_unbinned_observation()
    while obs != None:
        yield obs
        obs = self._next_unbinned_observation()
ctool._unbinned_observations   = _unbinned_observations
cscript._unbinned_observations = _unbinned_observations

# Run the script
def _run(self):
    if self._logDebug():
        self._log.cout(True)
    self._inc_running()
    self.process()
    self._dec_running()
csobservation.run = _run

# Execute the script
def _execute(self):
    self._inc_running()
    self._read_ahead(True)
    if self._logDebug():
        self._log.cout(True)
    self.process()
    self.save()
    self._dec_running()
csobservation.execute = _execute
%}
