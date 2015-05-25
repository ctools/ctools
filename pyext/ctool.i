/***************************************************************************
 *                        ctool - ctool base class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Juergen Knoedlseder                         *
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

    // Protected methods
    void                  init_members(void);
    void                  copy_members(const ctool& app);
    void                  free_members(void);
    const bool&           read_ahead(void) const;
    const GTimeReference& time_reference(void) const;
    GObservations         get_observations(const bool& get_response = true);

    // Protected methods that create objects from user parameters
    GEbounds        create_ebounds(void);
    GSkymap         create_map(const GObservations& obs);
    GCTAEventCube   create_cube(const GObservations& obs);
    GCTAObservation create_cta_obs(void);

    // Protected methods that check user parameters
    void            require_inobs(const std::string& method);

    // Protected support methods
    void            set_response(GObservations& obs);
    void            set_obs_response(GCTAObservation* obs);
    void            set_obs_bounds(GObservations& obs);
    GSkyDir         get_mean_pointing(const GObservations& obs);
    size_t          get_current_rss(void);
};


/***********************************************************************//**
 * @class cscript
 *
 * @brief cscript base class
 *
 * This is the base class from which all cscripts should derive. This
 * enabled using of ctool base class methods for generic parameter
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

