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
 * @brief Base class for ctools
 ***************************************************************************/
class ctool : public GApplication  {
public:
    // Constructors and destructors
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
    %rename(_read_ahead)               read_ahead;
    %rename(_time_reference)           time_reference;
    %rename(_get_observations)         get_observations;
    %rename(_setup_observations)       setup_observations;
    %rename(_setup_models)             setup_models;
    %rename(_create_ebounds)           create_ebounds;
    %rename(_create_map)               create_map;
    %rename(_create_cube)              create_cube;
    %rename(_create_cta_obs)           create_cta_obs;
    %rename(_require_inobs)            require_inobs;
    %rename(_require_inobs_nocube)     require_inobs_nocube;
    %rename(_get_roi)                  get_roi;
    %rename(_set_response)             set_response;
    %rename(_set_edisp)                set_edisp;
    %rename(_restore_edisp)            restore_edisp;
    %rename(_set_obs_response)         set_obs_response;
    %rename(_set_obs_bounds)           set_obs_bounds;
    %rename(_get_mean_pointing)        get_mean_pointing;
    %rename(_get_current_rss)          get_current_rss;
    %rename(_get_obs_header)           get_obs_header;
    %rename(_insert_energy_boundaries) insert_energy_boundaries;
    %rename(_cube_layer_usage)         cube_layer_usage;
    %rename(_is_valid_filename)        is_valid_filename;
    %rename(_warn_too_few_energies)    warn_too_few_energies;
    %rename(_warn_xml_suffix)          warn_xml_suffix;

    // Protected methods
    const bool&           read_ahead(void) const;
    const GTimeReference& time_reference(void) const;

    // Protected high-level setup methods
    void setup_observations(GObservations& obs, const bool& response = true,
                                                const bool& list = true,
                                                const bool& cube = true);
    void setup_models(GObservations& obs, const std::string& name = "");

    // Protected methods that create objects from user parameters
    GEbounds        create_ebounds(void);
    GSkyMap         create_map(const GObservations& obs);
    GCTAEventCube   create_cube(const GObservations& obs);
    GCTAObservation create_cta_obs(void);

    // Protected methods that check user parameters
    void            require_inobs(const std::string& method);
    void            require_inobs_nocube(const std::string& method);

    // Protected methods that extract user parameters
    GCTARoi         get_roi(void);

    // Protected support methods
    void              set_response(GObservations& obs);
    std::vector<bool> set_edisp(GObservations& obs, const bool& edisp) const;
    void              restore_edisp(GObservations& obs,
                                    const std::vector<bool>& edisp) const;
    void              set_obs_response(GCTAObservation* obs);
    void              set_obs_bounds(GObservations& obs);
    GObservations     get_observations(const bool& get_response = true);
    GSkyDir           get_mean_pointing(const GObservations& obs);
    size_t            get_current_rss(void);
    std::string       get_obs_header(const GObservation* obs) const;
    GEnergies         insert_energy_boundaries(const GEnergies&       energies,
                                               const GCTAObservation& obs);
    std::vector<bool> cube_layer_usage(const GEbounds& cube_ebounds,
                                       const GEbounds& list_ebounds) const;
    bool              is_valid_filename(const GFilename& filename) const;

    // Protected warning strings
    std::string     warn_too_few_energies(const GEnergies& energies) const;
    std::string     warn_xml_suffix(const GFilename& filename) const;
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
    void _log_observations(const int&           chatter,
                           const GObservations& obs,
                           const std::string&   what = "Observation") {
        self->log_observations(GChatter(chatter), obs, what);
    }
    void _log_models(const int&               chatter,
                           const GModels&     models,
                           const std::string& what = "Model") {
        self->log_models(GChatter(chatter), models, what);
    }
    void _read_ahead(const bool& flag) {
        self->m_read_ahead = flag;
        return;
    }
}


/***********************************************************************//**
 * @brief ctool and cscript base class Python extensions
 ***************************************************************************/
%pythoncode %{

# Initialise observation container from constructor arguments. In case that
# an observation container is provided as argument, this container will be
# used to initialise the class. The function returns a tuple containing the
# observation container and the eventually reduced argument list.
def _set_input_obs(self, argv):
    if len(argv) > 0 and isinstance(argv[0],gammalib.GObservations):
        obs  = argv[0]
        argv = argv[1:]
    else:      
        obs = gammalib.GObservations()
    return (obs, argv)
cscript._set_input_obs = _set_input_obs 

# Initialise application by calling the appropriate class constructor
def _init_cscript(self, argv):
    if len(argv) == 0:
        cscript.__init__(self, self._name, self._version)
    elif len(argv) == 1:
        cscript.__init__(self, self._name, self._version, *argv)
    else:
        raise TypeError('Invalid number of arguments given.')
    # Set logger properties
    self._log_header()
    self._log.date(True)
cscript._init_cscript = _init_cscript

# This function either returns or assigns user parameters in form of a
# dictionary. The assignment function takes a dictionary as argument and
# the function then loops over all keys in that dictionary and assigns
# the value to the parameter. For example
#
# >>> tool.pardict({'ra': 83.63, 'dec': 22.01})
#
# is equivalent to
#
# >>> tool['ra'] = 83.63
# >>> tool['dec'] = 22.01
#
# and assigns the 'ra' and 'dec' parameters of the tool. The return function
# takes no argument and returns a dictionary with the parameter names as
# keys and the current values as values. No parameter querying is performed.
# For example
#
# >>> d = tool.pardict()
#
# puts all parameters in the dictionary d.
def _pardict(self, *args):
    if len(args) == 0:
        d = {}
        for par in self._pars():
            if par.type() == 'b':
                v     = gammalib.tolower(par.current_value())
                value = (v == "yes" or v == "y" or v == "true" or v == "t")
                d[par.name()] = value
            elif par.type() == 'i':
                d[par.name()] = int(par.current_value())
            elif par.type() == 'r':
                d[par.name()] = float(par.current_value())
            else:
                d[par.name()] = str(par.current_value())
        return d
    elif len(args) == 1:
        for key in args[0]:
            self[key] = args[0][key]
    else:
        raise TypeError('pardict() takes 0 or 1 arguments (%d given)' % len(args))
ctool.pardict = _pardict
cscript.pardict = _pardict
%}
