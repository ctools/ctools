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

    // Initialise point source model
    init_testsource();

    // Get bins to be computed
    // To split the computation on several jobs,
    // the user is able to provide a subset of bins
    // which should be computed
    int binmin = (m_binmin == -1) ? 0 : m_binmin;
    int binmax = (m_binmax == -1) ? m_tsmap.npix() : m_binmax;

    GOptimizerLM* opt = (logTerse()) ? new GOptimizerLM(log)
        	                                     : new GOptimizerLM();

    // Get likelihood of null hypothesis if not given before
    if (m_logL0 == 0.0) {
    	m_obs.optimize(*opt);
    	m_logL0 = -(opt->value());
    }

    GModels models_orig = m_obs.models();
    GModels models = m_obs.models();

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

    	// append test source to model container
    	models.append(m_testsource);

    	// Assign models to observation
    	m_obs.models(models);

    	// Optimize observation container
    	m_obs.optimize(*opt);

    	// Retrieve the Likelihood value
    	double logL1 = -(opt->value());

    	// Compute TS value
    	double ts         = 2.0 * (logL1-m_logL0);

    	// Get test source model instance
    	GModels best_fit_model = m_obs.models();
    	GModel* testsource = best_fit_model["TestSource"];

    	// Assign values to the maps
    	m_tsmap(i) = ts;
    	m_fluxmap(i) = (*testsource)["Integral"].value();
    	m_indexmap(i) = (*testsource)["Index"].value();
    	m_statusmap(i) = 1.0;

    	// Remove test source from container
    	models.remove("TestSource");

    	// Bring models to initial state
    	m_obs.models(models_orig);

    } // endfor: looped over bins

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
    m_fluxmap.write(fitsfile);
    m_indexmap.write(fitsfile);

    // Change names of extensions
	fitsfile[0]->extname("TS MAP");
	fitsfile[1]->extname("FLUX MAP");
	fitsfile[2]->extname("SPECTRAL INDEX MAP");

	// Add computation log if not all bins are computed
    if (m_binmin != -1 || m_binmax != -1) {
    	m_statusmap.write(fitsfile);
    	fitsfile[3]->extname("STATUS MAP");
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
    m_fluxmap.clear();
    m_indexmap.clear();
    m_statusmap.clear();
    m_testsource.clear();
    m_source_integral = 1e-12;
    m_source_integral_min = 1e-16;
    m_source_integral_max = 1e-8;
    m_source_index = -2.5;
    m_source_index_min = - 5.0;
    m_source_index_max = -1.0;
    m_source_emin = GEnergy(1,"TeV");
    m_source_emax = GEnergy(100, "TeV");
    m_free_index = false;

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
    m_fluxmap       = app.m_fluxmap;
    m_indexmap       = app.m_indexmap;
    m_statusmap    = app.m_statusmap;
    m_testsource = app.m_testsource;
    m_source_integral = app.m_source_integral;
    m_source_integral_min = app.m_source_integral_min;
    m_source_integral_max = app.m_source_integral_max;
    m_source_index = app.m_source_integral_max;
    m_source_index_min = app.m_source_index_min;
    m_source_index_max = app.m_source_index_max;
    m_source_emin = app.m_source_emin;
    m_source_emax = app.m_source_emax;
    m_free_index = app.m_free_index;

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

    // Get map parameters
    m_xref = (*this)["xref"].real();
    m_yref = (*this)["yref"].real();
    m_proj     = (*this)["proj"].string();
    m_coordsys = (*this)["coordsys"].string();
    m_binsz    = (*this)["binsz"].real();
    m_nxpix    = (*this)["nxpix"].integer();
    m_nypix    = (*this)["nypix"].integer();

    // Get test source parameters
    m_source_integral = (*this)["source_integral"].real();
    m_source_integral_min = (*this)["source_integral_min"].real();
    m_source_integral_max = (*this)["source_integral_max"].real();
    m_source_index = (*this)["source_index"].real();
    m_source_index_min = (*this)["source_index_min"].real();
    m_source_index_max = (*this)["source_index_max"].real();
    m_source_emin = GEnergy((*this)["source_emin"].real(),"TeV");
    m_source_emax = GEnergy((*this)["source_emax"].real(),"TeV");
    m_free_index = (*this)["free_index"].boolean();

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
    m_fluxmap.clear();

    // Create skymap
    m_fluxmap = GSkymap(m_proj, m_coordsys,
                     m_xref, m_yref, -m_binsz, m_binsz,
                     m_nxpix, m_nypix, 1);

    // Initialise map information
    m_indexmap.clear();

    // Create skymap
    m_indexmap = GSkymap(m_proj, m_coordsys,
                     m_xref, m_yref, -m_binsz, m_binsz,
                     m_nxpix, m_nypix, 1);

    // Initialise map information
    m_statusmap.clear();

    // Create status map
    m_statusmap = GSkymap(m_proj, m_coordsys,
                     m_xref, m_yref, -m_binsz, m_binsz,
                     m_nxpix, m_nypix, 1);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise the test source
 *
 * Initialises the test point source model.
 ***************************************************************************/
void cttsmap::init_testsource(void)
{
	// Create spatial model
	GModelSpatialPointSource spat = GModelSpatialPointSource(m_tsmap.inx2dir(0));
	spat["RA"].fix();
	spat["DEC"].fix();

	// Create spectral model
	GModelSpectralPlaw2 plaw2 = GModelSpectralPlaw2(m_source_integral,m_source_index,m_source_emin,m_source_emax);

	// Create sky model
	m_testsource = GModelSky(spat,plaw2);
	m_testsource.name("TestSource");

	// Set source parameter boundaries
	m_testsource["Integral"].min(m_source_integral_min);
	m_testsource["Integral"].max(m_source_integral_max);
	m_testsource["Index"].min(m_source_index_min);
	m_testsource["Index"].max(m_source_index_max);

	// Free spectral index if required
	if (m_free_index) {
		m_testsource["Index"].free();
	}
	else {
		m_testsource["Index"].fix();
	}

    // Return
    return;
}

