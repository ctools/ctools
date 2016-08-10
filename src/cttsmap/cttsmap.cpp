/***************************************************************************
 *                    cttsmap - TS map calculation tool                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Michael Mayer                               *
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
 * @file cttsmap.hpp
 * @brief TS map calculation tool interface implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "cttsmap.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                          "cttsmap::get_parameters()"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
cttsmap::cttsmap(void) : ctool(CTTSMAP_NAME, CTTSMAP_VERSION)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observations constructor
 *
 * param[in] obs Observation container.
 *
 * This method creates an instance of the class by copying an existing
 * observations container.
 ***************************************************************************/
cttsmap::cttsmap(const GObservations& obs) :
         ctool(CTTSMAP_NAME, CTTSMAP_VERSION)
{
    // Initialise members
    init_members();

    // Set observations
    m_obs = obs;

    // Return
    return;
}



/***********************************************************************//**
 * @brief Command line constructor
 *
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 ***************************************************************************/
cttsmap::cttsmap(int argc, char *argv[]) :
         ctool(CTTSMAP_NAME, CTTSMAP_VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Application.
 ***************************************************************************/
cttsmap::cttsmap(const cttsmap& app) : ctool(app)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(app);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
cttsmap::~cttsmap(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] app Application.
 * @return Application.
 ***************************************************************************/
cttsmap& cttsmap::operator=(const cttsmap& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->ctool::operator=(app);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(app);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear cttsmap tool
 *
 * Clears cttsmap tool.
 ***************************************************************************/
void cttsmap::clear(void)
{
    // Free members
    free_members();
    this->ctool::free_members();

    // Clear base class (needed to conserve tool name and version)
    this->GApplication::clear();

    // Initialise members
    this->ctool::init_members();
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Computes the TS maps
 *
 * This method moves a source along a grid of coordinates and refits the
 * given models including this additional source. The fit parameters of the
 * source and its TS values are filled into respective position in the map.
 ***************************************************************************/
void cttsmap::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Set energy dispersion flags of all CTA observations and save old
    // values in save_edisp vector
    std::vector<bool> save_edisp = set_edisp(m_obs, m_apply_edisp);

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write header
    log_header1(TERSE, "Initialise TS map");

    // Get bins to be computed
    // To split the computation on several jobs,
    // the user is able to provide a subset of bins
    // which should be computed
    int binmin = (m_binmin == -1) ? 0 : m_binmin;
    int binmax = (m_binmax == -1) ? m_tsmap.npix() : m_binmax;

    // Initialise optimizer
    GOptimizerLM* opt = (logExplicit()) ? new GOptimizerLM(&log)
        	                            : new GOptimizerLM();

    // Store initial models
    GModels models_orig = m_obs.models();

    // Get model instance for further computations
    GModels models = m_obs.models();

    // Remove test source
    models.remove(m_srcname);

    // Get likelihood of null hypothesis if not given before
	if (m_logL0 == 0.0) {

        // Write header
        log_header1(TERSE, "Compute NULL Hypothesis for TS computation");

		// Compute likelihood without the test source
		m_obs.models(models);
		m_obs.optimize(*opt);
		m_logL0 = -(opt->value());

	}

    // Write header
    log_header1(TERSE, "Generate TS map");

    // Loop over grid positions
    for (int i = binmin; i < binmax; ++i) {

    	// Get the coordinate of current bin
    	GSkyDir bincentre = m_tsmap.inx2dir(i);

    	// Header for verbose logging
        log_header2(EXPLICIT, "Computing TS for bin number "+gammalib::str(i)+
                    " at "+bincentre.print());

    	// Add test source at current bin position
        if (m_testsource != NULL) {
            (*m_testsource)["RA"].value(bincentre.ra_deg());
            (*m_testsource)["DEC"].value(bincentre.dec_deg());
            models.append(*m_testsource);
        }

    	// Assign models to observations
    	m_obs.models(models);

    	// Optimize observation container
    	m_obs.optimize(*opt);

    	// Retrieve the Likelihood value
    	double logL1 = -(opt->value());

    	// Compute TS value
    	double ts = 2.0 * (logL1 - m_logL0);

    	// Log information
    	if (logExplicit()) {
    		log << " TS value ....: ";
            log << ts << std::endl;
    	}
        else if (logTerse()) {
    		log << "TS for bin number ";
            log << i;
            log << " at ";
            log << bincentre.print();
            log << ": ";
            log << ts << std::endl;
        }

    	// Get test source model instance
    	GModels best_fit_model = m_obs.models();
    	GModel* testsource     = best_fit_model[m_srcname];

    	// Assign values to the maps
    	m_tsmap(i) = ts;

    	// Extract fitted test source parameters
        if (m_testsource != NULL) {
            for (int j = 0; j < m_mapnames.size(); ++j) {
                m_maps[j](i) = (*testsource)[m_mapnames[j]].value();
            }
    	}

    	// Set status of bin to true
    	m_statusmap(i) = 1.0;

    	// Remove model from container
    	models.remove(m_srcname);

    } // endfor: looped over grid positions

    // Bring models to initial state
    m_obs.models(models_orig);

    // Restore energy dispersion flags of all CTA observations
    restore_edisp(m_obs, save_edisp);

    // Optionally publish TS map
    if (m_publish) {
        publish();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save maps
 *
 * This method saves the TS map and additional maps into a FITS file.
 ***************************************************************************/
void cttsmap::save(void)
{
    // Write header
    log_header1(TERSE, "Save TS map");

    // Get TS map filename
    m_outmap = (*this)["outmap"].filename();

    // Save only if filename is non-empty
    if (!m_outmap.is_empty() && !m_tsmap.is_empty()) {

        // Create fits file
        GFits fits;

        // Write the sky maps to the FITS file
        m_tsmap.write(fits);
        for (int i = 0; i < m_mapnames.size(); i++) {
            m_maps[i].write(fits);
        }

        // Set extension name for all maps
        for (int i = 0; i < m_mapnames.size(); i++) {
            fits[i+1]->extname(m_mapnames[i]);
        }

        // Add computation log if not all bins are computed
        if (m_binmin != -1 || m_binmax != -1) {
            m_statusmap.write(fits);
            fits[m_mapnames.size()+1]->extname("STATUS MAP");
        }

        // Save FITS file
        fits.saveto(m_outmap, clobber());

    } // endif: TS map filename was valid

    // Write into logger what has been done
    std::string fname = (m_outmap.is_empty()) ? "NONE" : m_outmap.url();
    if (m_tsmap.is_empty()) {
        fname.append(" (TS map is empty, no file created)");
    }
    log_value(NORMAL, "TS map file", fname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Publish TS map
 *
 * @param[in] name TS map name.
 ***************************************************************************/
void cttsmap::publish(const std::string& name)
{
    // Write header into logger
    log_header1(TERSE, "Publish TS map");

    // Set default name is user name is empty
    std::string user_name(name);
    if (user_name.empty()) {
        user_name = CTTSMAP_NAME;
    }

    // Write TS map name into logger
    log_value(NORMAL, "TS map name", user_name);

    // Publish TS map
    m_tsmap.publish(user_name);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void cttsmap::init_members(void)
{
    // Initialise members
    m_srcname.clear();
    m_outmap.clear();
    m_apply_edisp = false;
    m_publish     = false;

    // Initialise protected members
    m_obs.clear();
    m_binmin     = -1;
    m_binmax     = -1;
    m_logL0      = 0.0;
    m_tsmap.clear();
    m_statusmap.clear();
    m_mapnames.clear();
    m_maps.clear();
    m_testsource = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void cttsmap::copy_members(const cttsmap& app)
{
    // Copy attributes
    m_srcname     = app.m_srcname;
    m_outmap      = app.m_outmap;
    m_apply_edisp = app.m_apply_edisp;
    m_publish     = app.m_publish;

    // Copy protected members
    m_binmin    = app.m_binmin;
    m_binmax    = app.m_binmax;
    m_logL0     = app.m_logL0;
    m_obs       = app.m_obs;
    m_tsmap     = app.m_tsmap;
    m_mapnames  = app.m_mapnames;
    m_maps      = app.m_maps;
    m_statusmap = app.m_statusmap;

    // Clone protected members
    m_testsource = (app.m_testsource != NULL) ? app.m_testsource->clone() : NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void cttsmap::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * @exception GException::invalid_value
 *            Test source not found or no RA/DEC parameters found for test
 *            source.
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. Most parameters are only required if no observation exists so
 * far in the observation container. In this case, a single CTA observation
 * will be added to the container, using the definition provided in the
 * parameter file.
 ***************************************************************************/
void cttsmap::get_parameters(void)
{
    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {

        // Throw exception if no input observation file is given
        require_inobs(G_GET_PARAMETERS);

        // Get observation container
        m_obs = get_observations();

    } // endif: there was no observation in the container

    // ... otherwise add response information and energy boundaries in case
    // that they are missing
    else {
        setup_observations(m_obs);
    }

    // If there is are no models associated with the observations then
    // load now the model definition
    if (m_obs.models().size() == 0) {

        // Get models XML filename
        std::string filename = (*this)["inmodel"].filename();

        // Setup models for optimizing.
        m_obs.models(GModels(filename));

    } // endif: no models were associated with observations

    // Get name of test source and check container for this name
    m_srcname = (*this)["srcname"].string();
    if (!m_obs.models().contains(m_srcname)) {
        std::string msg = "Source \""+m_srcname+"\" not found in model "
                          "container. Please add a source with that name "
                          "or check for a possible typos.";
    	throw GException::invalid_value(G_GET_PARAMETERS, msg);
    }

    // Get Model and store it as protected member
    m_testsource = m_obs.models()[m_srcname]->clone();

    // Check if model has RA and DEC parameters which are necessary to
    // continue
    if (m_testsource != NULL) {
        if (!m_testsource->has_par("RA") || !m_testsource->has_par("DEC")) {
            std::string msg = "Source \""+m_srcname+"\" has no \"RA\" and "
                              "\"DEC\" parameters. Only sources with \"RA\" "
                              " and \"DEC\" parameters can be used as test "
                              "sources.";
            throw GException::invalid_value(G_GET_PARAMETERS, msg);
        }
    }

    // Fix the spatial parameters of the test source
    if (m_testsource != NULL) {
        (*m_testsource)["RA"].fix();
        (*m_testsource)["DEC"].fix();
    }

    // Create sky map based on task parameters
    GSkyMap map = create_map(m_obs);

    // Initialise maps from user parameters
    init_maps(map);

    // Read energy dispersion flag
    m_apply_edisp = (*this)["edisp"].boolean();

    // Get optional splitting parameters
    m_binmin  = (*this)["binmin"].integer();
    m_binmax  = (*this)["binmax"].integer();
    m_logL0   = (*this)["logL0"].real();
    m_publish = (*this)["publish"].boolean();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outmap = (*this)["outmap"].filename();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise skymaps
 *
 * Initialises skymaps that will contain map information.
 ***************************************************************************/
void cttsmap::init_maps(const GSkyMap& map)
{
    // Initialise map information
    m_tsmap.clear();

    // Create skymap
    m_tsmap = GSkyMap(map);

    // Initialise map information
	m_statusmap.clear();

	// Create status map
	m_statusmap = GSkyMap(map);

	// Initialise maps of free parameters
    if (m_testsource != NULL) {

        // Loop over all model parameters
        for (int i = 0; i < m_testsource->size(); i++) {
        
            // Skip fixed parameters
            if ((*m_testsource)[i].is_fixed()) {
                continue;
            }

            // Add sky map for free parameters. Note that push_back will
            // create a copy.
            m_maps.push_back(map);

            // Store parameter name
            m_mapnames.push_back((*m_testsource)[i].name());
        
        } // endfor: looped over all model parameters

	} // endif: test source was valid

    // Return
    return;
}
