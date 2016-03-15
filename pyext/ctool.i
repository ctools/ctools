/***************************************************************************
 *                        ctool - ctool base class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Juergen Knoedlseder                         *
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
 * @file ctool.hpp
 * @brief ctool base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#define SWIG
#include "ctool.hpp"


/***********************************************************************//**
 * @class cscript
 *
 * @brief cscript base class
 ***************************************************************************/
class cscript : public ctool  {
public:
    // Constructors and destructors
    cscript(void) : ctool() {}
    cscript(const std::string& name, const std::string& version) : ctool(name, version) {}
    cscript(const std::string& name, const std::string& version,
            int argc, char* argv[]) : ctool(name, version, argc, argv) {}
    cscript(const ctool& app) : ctool(app) {}
    virtual ~cscript(void) {}

    // Dummy methods (implementation makes this class non-abstract)
    virtual void clear(void) {
    }
    virtual void run(void) {
    }
    virtual void save(void) {
    }
};
%}

// Include (int ARGC, char **ARGV) typemap to allow passing command line
// arguments to GApplication constructor
%include "argcargv.i"


/***********************************************************************//**
 * @class ctool
 *
 * @brief ctool base class
 ***************************************************************************/
class ctool : public GApplication  {
public:
    // Constructors and destructors
    ctool(void);
    ctool(const std::string& name, const std::string& version);
    ctool(const std::string& name, const std::string& version,
          int ARGC, char **ARGV);
    ctool(const ctool& app);
    virtual ~ctool(void);


    // Pure virtual methods
    virtual void clear(void) = 0;
    virtual void run(void) = 0;
    virtual void save(void) = 0;

    // Public methods
    virtual void execute(void);

    // Make methods private in Python by prepending an underscore
    %rename(_read_ahead)           read_ahead() const;
    %rename(_time_reference)       time_reference() const;
    %rename(_get_observations)     get_observations(const bool& get_response = true);
    %rename(_create_ebounds)       create_ebounds();
    %rename(_create_map)           create_map(const GObservations& obs);
    %rename(_create_cube)          create_cube(const GObservations& obs);
    %rename(_create_cta_obs)       create_cta_obs();
    %rename(_require_inobs)        require_inobs(const std::string& method);
    %rename(_require_inobs_nolist) require_inobs_nolist(const std::string& method);
    %rename(_require_inobs_nocube) require_inobs_nocube(const std::string& method);
    %rename(_set_response)         set_response(GObservations& obs);
    %rename(_set_obs_response)     set_obs_response(GCTAObservation* obs);
    %rename(_set_obs_bounds)       set_obs_bounds(GObservations& obs);
    %rename(_get_mean_pointing)    get_mean_pointing(const GObservations& obs);
    %rename(_get_current_rss)      get_current_rss();
    %rename(_get_obs_header)       get_obs_header();

    // Protected methods
    const bool&           read_ahead(void) const;
    const GTimeReference& time_reference(void) const;
    GObservations         get_observations(const bool& get_response = true);

    // Protected methods that create objects from user parameters
    GEbounds        create_ebounds(void);
    GSkyMap         create_map(const GObservations& obs);
    GCTAEventCube   create_cube(const GObservations& obs);
    GCTAObservation create_cta_obs(void);

    // Protected methods that check user parameters
    void            require_inobs(const std::string& method);
    void            require_inobs_nolist(const std::string& method);
    void            require_inobs_nocube(const std::string& method);

    // Protected support methods
    void            set_response(GObservations& obs);
    void            set_obs_response(GCTAObservation* obs);
    void            set_obs_bounds(GObservations& obs);
    GSkyDir         get_mean_pointing(const GObservations& obs);
    size_t          get_current_rss(void);
    std::string     get_obs_header(const GObservation* obs);
};


/***********************************************************************//**
 * @class cscript
 *
 * @brief cscript base class
 *
 * This is the base class from which all cscripts should derive. This
 * enables using of ctool base class methods for generic parameter
 * handling.
 ***************************************************************************/
class cscript : public ctool  {
public:        
    // Constructors and destructors
    cscript(void);
    cscript(const std::string& name, const std::string& version);
    cscript(const std::string& name, const std::string& version,
            int ARGC, char **ARGV);
    cscript(const ctool& app);
    virtual ~cscript(void);

    // Methods
    virtual void clear(void);
    virtual void run(void);
    virtual void save(void);
};


/***********************************************************************//**
 * @brief ctool base class Python extensions
 ***************************************************************************/
%extend ctool {
    void _read_ahead(const bool& flag) {
        self->m_read_ahead = flag;
        return;
    }
}

/***********************************************************************//**
 * @brief cscript base class Python extensions
 ***************************************************************************/
%pythoncode %{

# Initialise observation container from constructor arguments. In case that
# an observation container is provided as argument, this container will be
# used to initialise the class.
def _set_input_obs(self, argv):
    if len(argv) > 0 and isinstance(argv[0],gammalib.GObservations):
        obs  = argv[0]
        argv = argv[1:]
    else:      
        obs = gammalib.GObservations()
    return obs
cscript._set_input_obs = _set_input_obs 

# Initialise application by calling the appropriate class constructor.
def _init_cscript(self, argv):
    if len(argv) == 0:
        cscript.__init__(self, self._name, self._version)
    elif len(argv) == 1:
        cscript.__init__(self, self._name, self._version, *argv)
    else:
        raise TypeError("Invalid number of arguments given.")
    # Set logger properties
    self._log_header()
    self._log.date(True)
cscript._init_cscript = _init_cscript

%}
