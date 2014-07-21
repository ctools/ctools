/***************************************************************************
 *                   ctbkgcube - CTA background cube tool                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Chia-Chun Lu                                     *
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
 * @file ctbkgcube.cpp
 * @brief CTA background cube tool implementation
 * @author Chia-Chun Lu
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctbkgcube.hpp"
#include "GTools.hpp"
#include "GWcs.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_EBOUNDS                              "ctbkgcube::get_ebounds()"
#define G_GET_PARAMETERS                        "ctbkgcube::get_parameters()"
#define G_SET_FROM_CNTMAP          "ctbkgcube::set_from_cntmap(std::string&)"

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
ctbkgcube::ctbkgcube(void) : GApplication(CTBKGCUBE_NAME, CTBKGCUBE_VERSION)
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
ctbkgcube::ctbkgcube(const GObservations& obs) :
           GApplication(CTBKGCUBE_NAME, CTBKGCUBE_VERSION)
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
ctbkgcube::ctbkgcube(int argc, char *argv[]) :
           GApplication(CTBKGCUBE_NAME, CTBKGCUBE_VERSION, argc, argv)
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
 * @param[in] app ctbkgcube application.
 ***************************************************************************/
ctbkgcube::ctbkgcube(const ctbkgcube& app) : GApplication(app)
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
ctbkgcube::~ctbkgcube(void)
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
 * @param[in] app ctbkgcube application.
 * @return Returns ctbkgcube application.
 ***************************************************************************/
ctbkgcube& ctbkgcube::operator=(const ctbkgcube& app)
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
void ctbkgcube::clear(void)
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
 * This is the main execution method of the ctbkgcube class. It is invoked
 * when the executable is called from command line. The method generates
 * the background cube and saves the result.
 ***************************************************************************/
void ctbkgcube::execute(void)
{
    // Signal that some parameters should be read ahead
    m_read_ahead = true;

    // Create the background cube
    run();

    // Save the background cube into the output FITS file
    save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generate the background cube(s).
 *
 * This method reads the task parameters from the parfile, sets up the
 * observation container, loops over all CTA observations in the container
 * and generates an background cube from the CTA observations.
 ***************************************************************************/
void ctbkgcube::run(void)
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

    // Set energy dispersion flag for all CTA observations and save old
    // values in save_edisp vector
    std::vector<bool> save_edisp;
    save_edisp.assign(m_obs.size(), false);
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            save_edisp[i] = obs->response()->apply_edisp();
            obs->response()->apply_edisp(m_apply_edisp);
        }
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
        log.header1("Generate background cube");
    }

    // If there are no models associated with the observations then
    // load now the model definition from the XML file
    if (m_obs.models().size() == 0) {
        m_bgmdl = (*this)["bgmdl"].filename();
        m_obs.models(GModels(m_bgmdl));
	
	// TODO: throw a message when no background definition file is given
    }

    // Loop over all observations in the container
    for (int i = 0; i < m_obs.size(); ++i) {
      
        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
 
	// Fill background cube
	fill_cube(obs, m_obs.models());

    }

    // Restore energy dispersion flag for all CTA observations
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            obs->response()->apply_edisp(save_edisp[i]);
        }
    }

    // Return
    return;
}

/***********************************************************************//**
 * @brief Generate background cube
 *
 * @param[in] obs CTA observation pointer.
 * @param[in] models Model container.
 *
 *
 ***************************************************************************/
void ctbkgcube::fill_cube(GCTAObservation* obs, const GModels& models)
{
  
    // Get pointing to do coordinate transformation
    GCTAPointing pnt = obs->pointing();
    GTime time = (*obs->events())[0]->time();
    double ontime = obs->ontime();
    // Continue only if observation pointer is valid
    if (obs != NULL) {

       // Initialise statistics
      double sum = 0.0;

      // Loop over all pixels in sky map
      for (int pixel = 0; pixel < m_bkgcube.npix(); ++pixel) {

	  GSkyDir dir     = m_bkgcube.inx2dir(pixel);
	  GCTAInstDir instdir = pnt.instdir(dir);
	  
	  // Loop over all background cube energy bins
	  for (int iebin = 0; iebin < m_ebounds.size(); ++iebin){
	 
            // Get event bin
            GCTAEventBin* bin = new GCTAEventBin();

            // Compute model value for event bin
            double model = 
                   models.eval(*(const_cast<const GCTAEventBin*>(bin)), *obs) *
                   bin->size();
	  
            // Store value
            bin->counts( bin->counts() + model );

            // Sum all events
            sum += model;
        }
      }

        // Log results
        if (logTerse()) {
            log << gammalib::parformat("Model events in cube");
            log << sum << std::endl;
        }

    } // endif: observation pointer was not valid

    // Return
    return;
}

/***********************************************************************//**
 * @brief Save background cube
 *
 ***************************************************************************/
void ctbkgcube::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save background cube");
    }

    // Get output filename
    m_outfile = (*this)["outfile"].filename();

    // Save background cube
    m_bkgcube.save(m_outfile, clobber());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. The parameters are read in the correct order.
 *
 * @todo Setup background cube from counts map
 ***************************************************************************/
void ctbkgcube::get_parameters(void)
{
    // If we do not have any observations in the container then get an
    // input file name or observation descriptor file
    if (m_obs.size() == 0) {
        get_obs();
    }

    // Make sure that response is set
    set_response();

    // Read energy dispersion flag
    m_apply_edisp = (*this)["edisp"].boolean();

    // If no counts map is specified then setup the background cube from
    // the user parameters
    std::string cntmap = (*this)["cntmap"].filename();
    if ((cntmap == "NONE") || (gammalib::strip_whitespace(cntmap) == "")) {
    
        // Get user parameters for counts map definition
        std::string wcs      = (*this)["proj"].string();
        std::string coordsys = (*this)["coordsys"].string();
        double      xref     = (*this)["xref"].real();
        double      yref     = (*this)["yref"].real();
        double      binsz    = (*this)["binsz"].real();
        int         nxpix    = (*this)["nxpix"].integer();
        int         nypix    = (*this)["nypix"].integer();

        // Get energy definition
        get_ebounds();

        // Define background cube
        m_bkgcube = GSkymap(wcs, coordsys, xref, yref,
			    -binsz, binsz, nxpix, nypix,
			    m_ebounds.size());
    }

    // ... otherwise setup the background cube from the counts map
    else {
    
        // Set background cube from counts map
        set_from_cntmap(cntmap);
    
    }

    // Read output filename (if needed)
    if (m_read_ahead) {
        m_outfile = (*this)["outfile"].filename();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get observation definition
 *
 * Get observation definition from the user parameters.
 ***************************************************************************/
void ctbkgcube::get_obs(void)
{
    // Get input filename
    std::string filename = (*this)["infile"].filename();

    // Try first to open as FITS file
    try {

        // Allocate CTA observation
        GCTAObservation obs;
        
        // Load input file in CTA observation
        obs.load(filename);

        // Append CTA observation to container
        m_obs.append(obs);

            
    }
        
    // ... otherwise try to open as XML file
    catch (GException::fits_open_error &e) {

        // Load observations from XML file
        m_obs.load(filename);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set observation response
 *
 * Set response for all observations that so have no response.
 ***************************************************************************/
void ctbkgcube::set_response(void)
{
    // Loop over all observations
    for (int i = 0; i < m_obs.size(); ++i) {

        // Is this observation a CTA observation?
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Yes ...
        if (obs != NULL) {

            // Set response if we don't have one
            if (!obs->hasresponse()) {

                // Load response information
                std::string database = (*this)["caldb"].string();
                std::string irf      = (*this)["irf"].string();

                // Set calibration database. If specified parameter is a
                // directory then use this as the pathname to the calibration
                // database. Otherwise interpret this as the instrument name,
                // the mission being "cta"
                GCaldb caldb;
                if (gammalib::dir_exists(database)) {
                    caldb.rootdir(database);
                }
                else {
                    caldb.open("cta", database);
                }

                // Set reponse
            	obs->response(irf, caldb);

            } // endif: observation already has a response

        } // endif: observation was a CTA observation

    } // endfor: looped over all observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get energy boundaries from parameters
 *
 * @exception GException::invalid_value
 *            Invalid extension name encountered.
 *
 * Get the energy boundaries from the user parameters.
 ***************************************************************************/
void ctbkgcube::get_ebounds(void)
{
    // Determine the energy binning alogrithm
    std::string ebinalg = (*this)["ebinalg"].string();

    // If we have the binning given by a file then try to get the boundaries
    // from that file
    if (ebinalg == "FILE") {
    
        // Get filename
        std::string filename = (*this)["ebinfile"].filename();

        // Open fits file to check which extension is given
        GFits file(filename);

        // Check first for EBOUNDS extension
        if (file.contains("EBOUNDS")) {
            file.close();
            m_ebounds.load(filename, "EBOUNDS");
        }

        // ... then check for ENERGYBINS extension
        else if (file.contains("ENERGYBINS")) {
            file.close();
            m_ebounds.load(filename, "ENERGYBINS");
        }

        // ... otherwise throw an exception
        else {
            file.close();
            std::string msg = "No extension with name \"EBOUNDS\" or"
                              " \"ENERGYBINS\" found in FITS file"
                              " \""+filename+"\".\n"
                              "An \"EBOUNDS\" or \"ENERGYBINS\" extension"
                              " is required if the parameter \"ebinalg\""
                              " is set to \"FILE\".";
            throw GException::invalid_value(G_GET_EBOUNDS, msg);
        }
    }
    
    // ... otherwise read emin, emax and nebins
    else {

        // Get the relevant parameters
    	double emin     = (*this)["emin"].real();
    	double emax     = (*this)["emax"].real();
    	int    enumbins = (*this)["enumbins"].integer();
        bool   log      = ((*this)["ebinalg"].string() == "LIN") ? false : true;

        // Create energy boundaries
        m_ebounds = GEbounds(enumbins,
                             GEnergy(emin, "TeV"),
                             GEnergy(emax, "TeV"),
                             log);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set background cube definition from counts map
 *
 * @exception GException::invalid_value
 *            Invalid counts map projection or invalid events encountered.
 *
 * Set background cube definition from counts map.
 ***************************************************************************/
void ctbkgcube::set_from_cntmap(const std::string& filename)
{
    // Allocate CTA observation
    GCTAObservation obs;
        
    // Load counts map in CTA observation
    obs.load(filename);

    // Set background cube from counts map
    const GCTAEventCube* cube = dynamic_cast<const GCTAEventCube*>(obs.events());

    // Continue only if cube is valid
    if (cube != NULL) {

        // Get sky map projection
        const GWcs* wcs = dynamic_cast<const GWcs*>(cube->map().projection());
        
        // Continue only if projection is valid
        if (wcs != NULL) {
            
            // Get user parameters for counts map definition
            std::string proj     = wcs->code();
            std::string coordsys = wcs->coordsys();
            double      xref     = wcs->crval(0);
            double      yref     = wcs->crval(1);
            double      dx       = wcs->cdelt(0);
            double      dy       = wcs->cdelt(1);
            int         nx       = cube->map().nx();
            int         ny       = cube->map().ny();

            // Get energy definition
            m_ebounds = cube->ebounds();

            // Define background cube
            m_bkgcube = GSkymap(proj, coordsys, xref, yref,
				dx, dy, nx, ny,
				m_ebounds.size());
        
        } // endif: WCS projection was valid

        // ... projection is not of WCS type
        else {
            std::string msg = "Counts map project is not of WCS type.";
            throw GException::invalid_value(G_SET_FROM_CNTMAP, msg);
        }

    } // endif: observation contained an events cube

    // ... there is not events cube
    else {
        std::string msg = "No events cube found in file \""
                          ""+filename+"\".";
        throw GException::invalid_value(G_SET_FROM_CNTMAP, msg);
    }

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
void ctbkgcube::init_members(void)
{
    // Initialise members
    m_outfile.clear();
    m_apply_edisp = false;

    // Initialise protected members
    m_read_ahead = false;
    m_obs.clear();
    m_bkgcube.clear();
    m_ebounds.clear();

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
void ctbkgcube::copy_members(const ctbkgcube& app)
{
    // Copy attributes
    m_outfile     = app.m_outfile;
    m_apply_edisp = app.m_apply_edisp;

    // Copy protected members
    m_read_ahead = app.m_read_ahead;
    m_obs        = app.m_obs;
    m_bkgcube    = app.m_bkgcube;
    m_ebounds    = app.m_ebounds;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctbkgcube::free_members(void)
{
    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}
