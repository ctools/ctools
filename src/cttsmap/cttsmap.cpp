/***************************************************************************
 *                      cttsmap - TS map computation tool                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file cttsmap.cpp
 * @brief TS map computation tool
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
#define G_GET_PARAMETERS                           "cttsmap::get_parameters()"

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
cttsmap::cttsmap(void) : GApplication(CTTSMAP_NAME, CTTSMAP_VERSION)
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

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
       GApplication(CTTSMAP_NAME, CTTSMAP_VERSION)
{
    // Initialise members
    init_members();

    // Set observations
    m_obs = obs;

    // Write header into logger
    log_header();

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
       GApplication(CTTSMAP_NAME, CTTSMAP_VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Application.
 ***************************************************************************/
cttsmap::cttsmap(const cttsmap& app) : GApplication(app)
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
 ***************************************************************************/
cttsmap& cttsmap::operator=(const cttsmap& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->GApplication::operator=(app);

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
 * @brief Clear instance
 ***************************************************************************/
void cttsmap::clear(void)
{
    // Free members
    free_members();
    this->GApplication::free_members();

    // Initialise members
    this->GApplication::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Execute application
 *
 * This is the main execution method of the cttsmap class. It is invoked when
 * the executable is called from command line.
 *
 * The method reads the task parameters, computed the TS values
 * and stores it into map(s), and writes the results into FITS files on disk.
 ***************************************************************************/
void cttsmap::execute(void)
{
    // Signal that some parameters should be read ahead
    m_read_ahead = true;

    // Compute the maps
    run();

    // Save the maps into FITS file
    save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Computes the TS maps
 *
 * This method moves a point-like source along a grid of coordinates
 * and refits the given models including this additional source.
 * The fit parameters of the source and its TS values are filled into
 * respective position in the map.
 ***************************************************************************/
void cttsmap::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Observations");
        }
        else {
            log.header1("Observation");
        }
        log << m_obs << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Bin observations");
        }
        else {
            log.header1("Bin observation");
        }
    }

    // Initialise maps
    init_maps();

    // Get bins to be computed
    // To split the computation on several jobs,
    // the user is able to provide a subset of bins
    // which should be computed
    int binmin = (m_binmin == -1) ? 0 : m_binmin;
    int binmax = (m_binmax == -1) ? m_tsmap.npix() : m_binmax;

    // Initialise optimizer
    GOptimizerLM* opt = (logTerse()) ? new GOptimizerLM(log)
        	                                     : new GOptimizerLM();

    // Store initial models
    GModels models_orig = m_obs.models();

    // Get model instance for further computations
    GModels models = m_obs.models();

    // Remove test source
    models.remove(m_srcname);

    // Get likelihood of null hypothesis if not given before
	if (m_logL0 == 0.0) {

		// Compute likelihood without the test source
		m_obs.models(models);
		m_obs.optimize(*opt);
		m_logL0 = -(opt->value());
	}

    // Loop over all observations in the container
    for (int i = binmin; i < binmax; ++i) {

    	// Get the coordinate of current bin
    	GSkyDir bincentre = m_tsmap.inx2dir(i);

    	// Verbosity option
    	if ((*this)["chatter"].integer()>2) {
    		std::cout<<"Computing TS for bin number "<<i<<" at "<<bincentre.print()<<std::endl;
    	}

    	// Set the point source at current bin position
    	m_testsource["RA"].value(bincentre.ra_deg());
    	m_testsource["DEC"].value(bincentre.dec_deg());

    	// Add test source to model container
    	models.append(m_testsource);

    	// Assign models to observations
    	m_obs.models(models);

    	// Optimize observation container
    	m_obs.optimize(*opt);

    	// Retrieve the Likelihood value
    	double logL1 = -(opt->value());

    	// Compute TS value
    	double ts         = 2.0 * (logL1-m_logL0);

    	// Get test source model instance
    	GModels best_fit_model = m_obs.models();
    	GModel* testsource = best_fit_model[m_srcname];

    	// Assign values to the maps
    	m_tsmap(i) = ts;

    	// Get fit values
    	for (int j = 0; j < m_mapnames.size(); j++) {

    		// get best fit values and assign to maps
    		m_maps[j](i) = (*testsource)[m_mapnames[j]].value();
    	}

    	// Set status of bin to true
    	m_statusmap(i) = 1.0;

    	// Remove model from container
    	models.remove(m_srcname);

    } // endfor: looped over bins

    // Bring models to initial state
    m_obs.models(models_orig);

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
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Save observations");
        }
        else {
            log.header1("Save observation");
        }
    }

    // Get output filename
    m_outfile = (*this)["outfile"].filename();

    // Create fits file instance
    GFits fitsfile = GFits();

    // Write the sky maps to the FITS file
    m_tsmap.write(fitsfile);

    for (int i = 0; i < m_mapnames.size(); i++) {
    	m_maps[i].write(fitsfile);
    }

    for (int i = 0; i < m_mapnames.size(); i++) {
    	fitsfile[i+1]->extname(m_mapnames[i]);
    }

	// Add computation log if not all bins are computed
    if (m_binmin != -1 || m_binmax != -1) {
    	m_statusmap.write(fitsfile);
    	fitsfile[m_mapnames.size()+1]->extname("STATUS MAP");
    }

    // Save FITS file
    fitsfile.saveto(m_outfile,true);

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
    m_evfile.clear();
    m_outfile.clear();
    m_proj.clear();
    m_coordsys.clear();
    m_xref     = 0.0;
    m_yref     = 0.0;
    m_binsz    = 0.0;
    m_nxpix    = 0;
    m_nypix    = 0;

    // Initialise protected members
    m_obs.clear();
    m_read_ahead = false;
    m_binmin = -1;
    m_binmax = -1;
    m_logL0 = 0.0;
    m_tsmap.clear();
    m_maps.clear();
    m_mapnames.clear();
    m_statusmap.clear();
    m_testsource.clear();

    // Set logger properties
    log.date(true);

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
    m_evfile   = app.m_evfile;
    m_outfile  = app.m_outfile;
    m_proj     = app.m_proj;
    m_coordsys = app.m_coordsys;
    m_xref     = app.m_xref;
    m_yref     = app.m_yref;
    m_binsz    = app.m_binsz;
    m_nxpix    = app.m_nxpix;
    m_nypix    = app.m_nypix;

    // Copy protected members
    m_binmin = app.m_binmin;
    m_binmax = app.m_binmax;
    m_logL0 = app.m_logL0;
    m_obs        = app.m_obs;
    m_read_ahead = app.m_read_ahead;
    m_tsmap       = app.m_tsmap;
    m_mapnames       = app.m_mapnames;
    m_maps       = app.m_maps;
    m_statusmap    = app.m_statusmap;
    m_testsource    = app.m_testsource;
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void cttsmap::free_members(void)
{
    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. Most parameters are only required if no observation exists so
 * far in the observation container. In this case, a single CTA observation
 * will be added to the container, using the definition provided in the
 * parameter file.
 ***************************************************************************/
void cttsmap::get_parameters(void)
{
    // If there are no observations in container then add a single CTA
    // observation using the parameters from the parameter file
    if (m_obs.size() == 0) {

        // Get name of CTA events file
        m_evfile = (*this)["evfile"].filename();

        // Allocate CTA observation
        GCTAObservation obs;

        // Try first to open as FITS file
        try {

            // Load event list in CTA observation
            obs.load(m_evfile);

            // Append CTA observation to container
            m_obs.append(obs);

        }

        // ... otherwise try to open as XML file
        catch (GException::fits_open_error &e) {

            // Load observations from XML file
            m_obs.load(m_evfile);

        }

        m_obs.models((*this)["srcmdl"].filename());

    } // endif: there was no observation in the container

    // Get name of test source and check container for this name
    m_srcname = (*this)["srcname"].string();
    if (!m_obs.models().contains(m_srcname)) {
    	throw GException::model_not_found(G_GET_PARAMETERS, m_srcname);
    }

    // Get Model and store it as protected member
    GModel* model = m_obs.models()[m_srcname]->clone();
    m_testsource = *dynamic_cast<GModelSky*>(model);

    // Check if model has RA and DEC parameters which are necessary to continue
    if (!m_testsource.has_par("RA") || !m_testsource.has_par("DEC")) {

    	// Is this the correct exception? Do we need a new one?
    	throw GException::model_invalid(G_GET_PARAMETERS, m_srcname, "RA and DEC parameters required! ");
    }

    // Fix the spatial parameters
    m_testsource["RA"].fix();
    m_testsource["DEC"].fix();

    // Get map parameters
    m_xref = (*this)["xref"].real();
    m_yref = (*this)["yref"].real();
    m_proj     = (*this)["proj"].string();
    m_coordsys = (*this)["coordsys"].string();
    m_binsz    = (*this)["binsz"].real();
    m_nxpix    = (*this)["nxpix"].integer();
    m_nypix    = (*this)["nypix"].integer();

    // Get optional splitting parameters
    m_binmin = (*this)["binmin"].integer();
    m_binmax = (*this)["binmax"].integer();
    m_logL0 = (*this)["logL0"].real();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (m_read_ahead) {
        m_outfile = (*this)["outfile"].filename();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise map information
 *
 * Initialises the skymaps.
 ***************************************************************************/
void cttsmap::init_maps(void)
{
    // Initialise map information
    m_tsmap.clear();

    // Create skymap
    m_tsmap = GSkymap(m_proj, m_coordsys,
                     m_xref, m_yref, -m_binsz, m_binsz,
                     m_nxpix, m_nypix, 1);

    // Initialise map information
	m_statusmap.clear();

	// Create status map
	m_statusmap = GSkymap(m_proj, m_coordsys,
					 m_xref, m_yref, -m_binsz, m_binsz,
					 m_nxpix, m_nypix, 1);

	// Initialise maps of other free parameters
	for (int i = 0; i < m_testsource.size(); i++) {
		if (m_testsource[i].is_fixed()) {
			continue;
		}

		// create sky map
		GSkymap map = GSkymap(m_proj, m_coordsys,
				 m_xref, m_yref, -m_binsz, m_binsz,
				 m_nxpix, m_nypix, 1);
		m_maps.push_back(map);

		// Store name of parameter
		m_mapnames.push_back(m_testsource[i].name());
	}

    // Return
    return;
}

