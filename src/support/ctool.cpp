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
 * @brief ctool base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>         // std::getenv() function
#include <cstdio>          // std::fopen(), etc. functions
#include "ctool.hpp"
#include "GTools.hpp"

/* __ Includes for memory usage determination ____________________________ */
#if defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>
#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>
#endif
#endif


/* __ Method name definitions ____________________________________________ */
#define G_GET_MEAN_POINTING        "ctool::get_mean_pointing(GObservations&)"
#define G_CREATE_EBOUNDS                            "ctool::create_ebounds()"
#define G_PROVIDE_HELP                                "ctool::provide_help()"

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
ctool::ctool(void) : GApplication()
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Name constructor
 *
 * @param[in] name Application name.
 * @param[in] version Application version.
 ***************************************************************************/
ctool::ctool(const std::string& name, const std::string& version) :
       GApplication(name, version)
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Command line constructor
 *
 * @param[in] name Application name.
 * @param[in] version Application version.
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 ***************************************************************************/
ctool::ctool(const std::string& name, const std::string& version,
             int argc, char *argv[]) : 
       GApplication(name, version, argc, argv)
{
    // Catch --help option before doing anything else
    if (need_help()) {
        provide_help();
        exit(0);
    }
    
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
ctool::ctool(const ctool& app) : GApplication(app)
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
ctool::~ctool(void)
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
ctool& ctool::operator=(const ctool& app)
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
 * @brief Execute application
 *
 * This is the main execution method of a ctool. The method is invoked when
 * the executable is called from command line.
 ***************************************************************************/
void ctool::execute(void)
{
    // Signal that some parameters should be read ahead
    m_read_ahead = true;

    // Run the tool
    run();

    // Save the results
    save();

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
void ctool::init_members(void)
{
    // Initialise members
    m_read_ahead = false;
    m_use_xml    = false;

    // Set CTA time reference. G_CTA_MJDREF is the CTA reference MJD,
    // which is defined in GCTALib.hpp. This is somehow a kluge. We need
    // a better mechanism to implement the CTA reference MJD.
    m_cta_ref.set(G_CTA_MJDREF, "s", "TT", "LOCAL");

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
void ctool::copy_members(const ctool& app)
{
    // Copy members
    m_read_ahead = app.m_read_ahead;
    m_use_xml    = app.m_use_xml;
    m_cta_ref    = app.m_cta_ref;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctool::free_members(void)
{
    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get observation container
 *
 * @param[in] get_response Indicates whether response information should
 *                         be loaded (default: true)
 *
 * Get an observation container according to the user parameters. The method
 * supports loading of a individual FITS file or an observation definition
 * file in XML format. If the input filename is empty, parameters are read
 * to build a CTA observation from scratch.
 ***************************************************************************/
GObservations ctool::get_observations(const bool& get_response)
{
    // Initialise empty observation container
    GObservations obs;

    // Get the filename from the input parameters
    std::string filename = (*this)["inobs"].filename();

    // If no observation definition file has been specified then read all
    // parameters that are necessary to create an observation from scratch
    if ((filename == "NONE") || (gammalib::strip_whitespace(filename) == "")) {

        // Setup a new CTA observation
        GCTAObservation cta_obs = create_cta_obs();

        // Get response if required
        if (get_response) {

            // Set response
            set_obs_response(&cta_obs);

        } // endif: response was required

       // Append observation to container
       obs.append(cta_obs);

    } // endif: filename was "NONE" or ""

    // ... otherwise we have a file name
    else {

        // If file is a FITS file then create an empty CTA observation
        // and load file into observation
        if (GFilename(filename).is_fits()) {

            // Allocate empty CTA observation
            GCTAObservation cta_obs;

            // Load data
            cta_obs.load(filename);

            // Get response if required
            if (get_response) {

                // Set response
                set_obs_response(&cta_obs);

            } // endif: response was required

            // Append observation to container
            obs.append(cta_obs);

            // Signal that no XML file should be used for storage
            m_use_xml = false;

        }

        // ... otherwise load file into observation container
        else {

            // Load observations from XML file
            obs.load(filename);

            // Get response if required
            if (get_response) {

                // For all observations that have no response, set the response
                // from the task parameters
                set_response(obs);

            } // endif: response was required

            // Set observation boundary parameters (emin, emax, rad)
            set_obs_bounds(obs);

            // Signal that XML file should be used for storage
            m_use_xml = true;

        } // endelse: file was an XML file

    }

    // Return observation container
    return obs;

}


/***********************************************************************//**
 * @brief Create energy boundaries from user parameters
 *
 * @exception GException::invalid_value
 *            No valid energy boundary extension found.
 *
 * Get the energy boundaries according to the user parameters. The method
 * supports loading of energy boundary information from the EBOUNDS or
 * ENERGYBINS extension, or setting energy boundaries using a linear
 * or logarithmical spacing.
 *
 * The following parameters are read:
 *      ebinalg - Energy binning algorithm
 *      emin - Minimum energy (if ebinalg != FILE)
 *      emax - Maximum energy (if ebinalg != FILE)
 *      enumbins - Number of energy bins (if ebinalg != FILE)
 ***************************************************************************/
GEbounds ctool::create_ebounds(void)
{
    // Allocate energy boundaries
    GEbounds ebounds;

    // Get energy binning algorithm
    std::string ebinalg = (*this)["ebinalg"].string();

    // If energy binning algorithm is of type "FILE" (case sensitive), then
    // read energy boundaries from FITS file ...
    if (ebinalg == "FILE") {

        // Get filename
        std::string ebinfile = (*this)["ebinfile"].filename();

        // Create filename instance
        GFilename fname(ebinfile);

        // Check if extension name was provided
        if (!fname.has_extname()) {

            // Open energy boundary file using the EBOUNDS or ENERGYBINS
            // extension. Throw an exception if opening fails.
            GFits file(ebinfile);
            if (file.contains("EBOUNDS")) {
                file.close();
                ebounds.load(ebinfile);
            }
            else if (file.contains("ENERGYBINS")) {
                file.close();
                ebinfile += "[ENERGYBINS]";
                ebounds.load(ebinfile);
            }
            else {
                file.close();
                std::string msg = "No extension with name \"EBOUNDS\" or"
                                  " \"ENERGYBINS\" found in FITS file"
                                  " \""+ebinfile+"\".\n"
                                  "An \"EBOUNDS\" or \"ENERGYBINS\" extension"
                                  " is required if the parameter \"ebinalg\""
                                  " is set to \"FILE\".";
                throw GException::invalid_value(G_CREATE_EBOUNDS, msg);
            }
        }
        else {
            // Load ebounds from filename including extension
            ebounds.load(ebinfile);
        }
        
    } // endif: ebinalg was "FILE"

    // ... otherwise use a linear or a logarithmically-spaced energy binning
    else {

        // Get task parameters
    	double emin     = (*this)["emin"].real();
    	double emax     = (*this)["emax"].real();
    	int    enumbins = (*this)["enumbins"].integer();

        // Initialise log mode for ebinning
        bool log = true;

        // Check if algorithm is linear
        if (ebinalg == "LIN") {
            log = false;
        }

        // Setup energy bins
        ebounds = GEbounds(enumbins, GEnergy(emin, "TeV"),
                                     GEnergy(emax, "TeV"), log);

    } // endelse: ebinalg was not "FILE"

    // Return energy boundaries
    return ebounds;
}


/***********************************************************************//**
 * @brief Create a skymap from user parameters
 *
 * @param[in] obs Observation container.
 *
 * Creates an empty skymap from input user parameters. The reference
 * coordinate is extracted from the mean pointing of the CTA observations
 * in the observation container (see get_mean_pointing()).
 *
 * The following parameters are read:
 *      usepnt - Use pointing for reference coordinate?
 *      xref - Reference in Right Ascension or longitude (if usepnt=no)
 *      yref - Reference in Declination or latitude (if usepnt=no)
 *      proj - Coordinate projection
 *      coordsys - Coordinate system
 *      binsz - Bin size (deg)
 *      nxpix - Number of pixels in Right Ascension or longitude
 *      nypix - Number of pixels in Declination or latitude
 ***************************************************************************/
GSkyMap ctool::create_map(const GObservations& obs)
{
    // Read task parameters
    double xref   = 0.0;
    double yref   = 0.0;
    bool   usepnt = (*this)["usepnt"].boolean();
    if (!usepnt) {
        xref = (*this)["xref"].real();
        yref = (*this)["yref"].real();
    }
    std::string proj     = (*this)["proj"].string();
    std::string coordsys = (*this)["coordsys"].string();
    double      binsz    = (*this)["binsz"].real();
    int         nxpix    = (*this)["nxpix"].integer();
    int         nypix    = (*this)["nypix"].integer();

    // If requested, get pointing from observations
    if (usepnt) {
    
        // Get mean pointing
        GSkyDir pnt = get_mean_pointing(obs);

        // Use correct coordinate system
        if (gammalib::toupper(coordsys) == "GAL") {
            xref = pnt.l_deg();
            yref = pnt.b_deg();
        }
        else {
            xref = pnt.ra_deg();
            yref = pnt.dec_deg();
        }

    } // endif: got mean pointing as reference

    // Initialise sky map
    GSkyMap map = GSkyMap(proj, coordsys, xref, yref, -binsz, binsz, 
                          nxpix, nypix, 1);

    // Return sky map
    return map;
}


/***********************************************************************//**
 * @brief Create a CTA event cube from user parameters
 *
 * @param[in] obs Observation container.
 *
 * Creates an empty CTA event cube from input user parameters. The reference
 * coordinate is extracted from the mean pointing of the CTA observations
 * in the observation container (see get_mean_pointing()). This method will
 * be used if a ctool has the input parameter usepnt=true. The method
 * appends a dummy GTI to the event cube.
 *
 * The following parameters are read:
 *      proj - Coordinate projection
 *      coordsys - Coordinate system
 *      binsz - Bin size (deg)
 *      nxpix - Number of pixels in Right Ascension or longitude
 *      nypix - Number of pixels in Declination or latitude
 *      ebinalg - Energy binning algorithm
 *      emin - Minimum energy (if ebinalg != FILE)
 *      emax - Maximum energy (if ebinalg != FILE)
 *      enumbins - Number of energy bins (if ebinalg != FILE)
 ***************************************************************************/
GCTAEventCube ctool::create_cube(const GObservations& obs)
{
    // Get skymap
    GSkyMap map = create_map(obs);

    // Set energy boundaries
    GEbounds ebounds = create_ebounds();

    // Set dummy GTI that is needed for event cube creation
    GGti gti;
    gti.append(GTime(0.0), GTime(0.1234));

    // Extend skymap to the requested number of maps
    map.nmaps(ebounds.size());

    // Initialise event cube from all parameters
    GCTAEventCube cube = GCTAEventCube(map, ebounds, gti);

    // Return event cube
    return cube;
}


/***********************************************************************//**
 * @brief Create a CTA observation from user parameters
 *
 * Creates an empty CTA observation from user parameters. An empty event list
 * including RoI, GTI and energy boundary information is attached to the
 * observation. The method also sets the pointing direction using the ra and
 * dec parameter, the ROI based on ra, dec and rad, a single GTI based on
 * tmin and tmax, and a single energy boundary based on emin and emax. The
 * method furthermore sets the ontime, livetime and deadtime correction
 * factor.
 *
 * The following parameters are read:
 *      ra - Right Ascension of pointing and RoI centre (deg)
 *      dec - Declination of pointing and RoI centre (deg)
 *      rad - Radius of RoI (deg)
 *      deadc - Deadtime correction factor
 *      tmin - Start time
 *      tmax - Stop time
 *      emin - Minimum energy (TeV)
 *      emax - Maximum energy (TeV)
 ***************************************************************************/
GCTAObservation ctool::create_cta_obs(void)
{
    // Get CTA observation parameters
    double ra    = (*this)["ra"].real();
    double dec   = (*this)["dec"].real();
    double rad   = (*this)["rad"].real();
    double deadc = (*this)["deadc"].real();
    double t_min = (*this)["tmin"].real();
    double t_max = (*this)["tmax"].real();
    double e_min = (*this)["emin"].real();
    double e_max = (*this)["emax"].real();

    // Allocate CTA observation and empty event list
    GCTAObservation obs;
    GCTAEventList   list;

    // Set pointing direction
    GCTAPointing pnt;
    GSkyDir      skydir;
    skydir.radec_deg(ra, dec);
    pnt.dir(skydir);

    // Set ROI
    GCTAInstDir instdir(skydir);
    GCTARoi     roi(instdir, rad);

    // Set GTI
    GGti  gti(m_cta_ref);
    GTime tstart;
    GTime tstop;
    tstart.set(t_min, m_cta_ref);
    tstop.set(t_max, m_cta_ref);
    gti.append(tstart, tstop);

    // Set energy boundaries
    GEbounds ebounds;
    GEnergy  emin;
    GEnergy  emax;
    emin.TeV(e_min);
    emax.TeV(e_max);
    ebounds.append(emin, emax);

    // Set CTA event list attributes
    list.roi(roi);
    list.gti(gti);
    list.ebounds(ebounds);

    // Attach empty event list to CTA observation
    obs.events(list);

    // Set observation ontime, livetime and deadtime correction factor
    obs.pointing(pnt);
    obs.ontime(gti.ontime());
    obs.livetime(gti.ontime()*deadc);
    obs.deadc(deadc);

    // Return CTA observation
    return obs;

}


/***********************************************************************//**
 * @brief Throws exception if inobs parameter is not valid
 *
 * @param[in] method Method name.
 *
 * Throw an exception if the inobs parameter is either NONE or an empty
 * string.
 ***************************************************************************/
void ctool::require_inobs(const std::string& method)
{
    // Get inobs filename
    std::string filename = (*this)["inobs"].filename();

    // Throw exception if no infile is given
    if (filename == "NONE" || gammalib::strip_whitespace(filename) == "") {
        std::string msg = "A valid file needs to be specified for the "
                          "\"inobs\" parameter, yet \""+filename+
                          "\" was given."
                          " Specify a valid observation definition or "
                          "FITS file to proceed.";
        throw GException::invalid_value(method, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Throws exception if inobs parameter is an event list
 *
 * @param[in] method Method name.
 *
 * Throw an exception if the inobs parameter is an event list.
 ***************************************************************************/
void ctool::require_inobs_nolist(const std::string& method)
{
    // Get inobs filename
    GFilename filename = (*this)["inobs"].filename();

    // Continue only if we have a FITS file
    if (filename.is_fits()) {

        // Signal no list
        bool is_list = false;
    
        // Try loading file as counts cube. If this is successful then
        // throw an exception
        try {

            // Load list from file
            GCTAEventList list(filename);

            // If we're still alive then signal that we have a list
            is_list = true;

        }

        // Catch any exceptions
        catch (...) {
            ;
        }

        // If we have an event list then throw an exception
        if (is_list) {
            std::string msg = "An event list has been specified for the "
                              "\"inobs\" parameter, yet no event list "
                              "can be specified as input observation."
                              " Instead, specify a counts cube or an "
                              "observation definition file.";
            throw GException::invalid_value(method, msg);
        }
    
    } // endif: we had a FITS file
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Throws exception if inobs parameter is a counts cube
 *
 * @param[in] method Method name.
 *
 * Throw an exception if the inobs parameter is a counts cube.
 ***************************************************************************/
void ctool::require_inobs_nocube(const std::string& method)
{
    // Get inobs filename
    GFilename filename = (*this)["inobs"].filename();

    // Continue only if we have a FITS file
    if (filename.is_fits()) {

        // Signal no cube
        bool is_cube = false;
    
        // Try loading file as counts cube. If this is successful then
        // throw an exception
        try {

            // Load cube from file
            GCTAEventCube cube(filename);

            // If we're still alive then signal that we have a cube
            is_cube = true;

        }

        // Catch any exceptions
        catch (...) {
            ;
        }

        // If we have a counts cube then throw an exception
        if (is_cube) {
            std::string msg = "A counts cube has been specified for the "
                              "\"inobs\" parameter, yet no counts cube "
                              "can be specified as input observation."
                              " Instead, specify an event list or an "
                              "observation definition file.";
            throw GException::invalid_value(method, msg);
        }
    
    } // endif: we had a FITS file
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set response for all CTA observations in container
 *
 * @param[in,out] obs Observation container
 *
 * Set the response for a CTA observations in the container that so far have
 * no response using the "caldb" and "irf" task parameters.
 *
 * The following parameters are read:
 *      caldb - Calibration database
 *      irf - Instrument response function
 ***************************************************************************/
void ctool::set_response(GObservations& obs)
{
    // Setup response for all observations
    for (int i = 0; i < obs.size(); ++i) {

        // Is this observation a CTA observation?
        GCTAObservation* cta = dynamic_cast<GCTAObservation*>(obs[i]);

        // Yes ...
        if (cta != NULL) {

            // Set response if we don't have one
            if (!cta->has_response()) {
                set_obs_response(cta);
            }

        } // endif: observation was a CTA observation

    } // endfor: looped over all observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set response for CTA observation
 *
 * @param[in,out] obs CTA observation
 *
 * Set the response for a CTA observation. If the CTA observation contains
 * a counts cube the method first attempts to set the response information
 * using the @a expcube and @a psfcube user parameters. If this is not
 * successful (for example because the parameters do not exist or their
 * value is NONE or blank), the method will set the response information
 * from the @a caldb and @a irf user parameters.
 *
 * The following parameters are read:
 *      expcube - Exposure cube file (optional)
 *      psfcube - PSF cube file (optional)
 *      caldb - Calibration database (optional)
 *      irf - Instrument response function (optional)
 ***************************************************************************/
void ctool::set_obs_response(GCTAObservation* obs)
{
    // Initialise response flag
    bool has_response = false;

    // If the observation contains a counts cube, then first check whether
    // the expcube and psfcube parameters exist and are not NONE
    if (dynamic_cast<const GCTAEventCube*>(obs->events()) != NULL) {
        if (has_par("expcube")   && has_par("psfcube") &&
            has_par("edispcube") && has_par("bkgcube")) {

            // Get filenames
            std::string expcube   = (*this)["expcube"].filename();
            std::string psfcube   = (*this)["psfcube"].filename();
            std::string edispcube = (*this)["edispcube"].filename();
            std::string bkgcube   = (*this)["bkgcube"].filename();

            // Extract response information if available
            if ((expcube != "NONE") && (psfcube != "NONE") &&
                (bkgcube != "NONE") &&
                (gammalib::strip_whitespace(expcube) != "") &&
                (gammalib::strip_whitespace(psfcube) != "") &&
                (gammalib::strip_whitespace(bkgcube) != "")) {

                // Get exposure, PSF and background cubes
                GCTACubeExposure   exposure(expcube);
                GCTACubePsf        psf(psfcube);
                GCTACubeBackground background(bkgcube);

                // Optionally load energy dispersion cube
                if ((edispcube != "NONE") &&
                    (gammalib::strip_whitespace(edispcube) != "")) {
                	GCTACubeEdisp edisp(edispcube);
                	obs->response(exposure, psf, edisp, background);
                }
                else {
					// Set reponse
					obs->response(exposure, psf, background);
                }

                // Signal response availability
                has_response = true;

            } // endif: filenames were available

        } // endif: expcube and psfcube parameters exist
    } // endif: observation contains a counts cube

    // If we have not yet response information then get it now
    // from the caldb and irf parameters
    if (!has_response) {

        // Load response information
        std::string database = (*this)["caldb"].string();
        std::string irf      = (*this)["irf"].string();

        // Create an XML element containing the database and IRF name. This
        // kluge will make sure that the information is later written
        // to the observation definition XML file, in case an observation
        // definition XML file is written.
        std::string parameter = "parameter name=\"Calibration\""
                                " database=\""+database+"\""
                                " response=\""+irf+"\"";

        GXmlElement xml;
        xml.append(parameter);

        // Create CTA response
        GCTAResponseIrf response(xml);

        // Attach response to observation
        obs->response(response);
        
        // Signal response availability
        has_response = true;

    } // endif: no response information was available

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set observation boundaries
 *
 * @param[in,out] obs Observation container
 ***************************************************************************/
void ctool::set_obs_bounds(GObservations& obs)
{
    // Setup response for all observations
    for (int i = 0; i < obs.size(); ++i) {

        // Is this observation a CTA observation?
        GCTAObservation* cta = dynamic_cast<GCTAObservation*>(obs[i]);

        // Continue only if observation is CTA and has events
        if ((cta != NULL) && (cta->has_events())) {

            // Get pointer on event list
            GCTAEventList* list = const_cast<GCTAEventList*>(dynamic_cast<const GCTAEventList*>(cta->events()));

            // Continue only if it's valid
            if (list != NULL) {

                // If there are no energy boundaries then read the user
                // parameters and add them
                if (list->ebounds().is_empty()) {

                    // If there are "emin" and "emax" parameters then use
                    // them
                    if (has_par("emin") && has_par("emax")) {
                        double emin = (*this)["emin"].real();
                        double emax = (*this)["emax"].real();
                        GEbounds ebounds(GEnergy(emin, "TeV"),
                                         GEnergy(emax, "TeV"));
                        list->ebounds(ebounds);
                    }
                    
                }

                // If there is no ROI then read the user parameters and add
                // them
                if (list->roi().radius() == 0) {
                    
                    // If there is a "rad" parameter then use it
                    if (has_par("rad")) {
                        double rad = (*this)["rad"].real();
                        GCTARoi roi(GCTAInstDir(cta->pointing().dir()), rad);
                        list->roi(roi);
                    }

                }

            } // endif: list was valid

        } // endif: observation was CTA
        
    } // endfor: looped over observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Derives mean pointing from CTA observations
 *
 * @param[in] obs Observation container.
 * @return Pointing direction.
 *
 * Computes the mean pointing direction from all CTA observations in the
 * observation container.
 *
 * @todo This method does not work properly for wrap arounds in Right
 *       Ascension. It certainly will also lead to bad results near the
 *       celestial pole.
 ***************************************************************************/
GSkyDir ctool::get_mean_pointing(const GObservations& obs)
{
    // Initialise values
    double xref  = 0.0;
    double yref  = 0.0;
    int    n_pnt = 0;

    // Loop over all observations to compute mean pointings
    for (int i = 0; i < obs.size(); ++i) {

        // If we have a CTA observation then extract the pointing direction
        const GCTAObservation* cta_obs = dynamic_cast<const GCTAObservation*>(obs[i]);
        if (cta_obs != NULL) {

           // Get pointing
           const GCTAPointing& pnt = cta_obs->pointing();

           // Add to coordinates
           xref  += pnt.dir().ra_deg();
           yref  += pnt.dir().dec_deg();
           n_pnt += 1;

       } // endif: cta observation was valid

    } // endfor: loop over all observations

    // Signal if no pointing is found
    if (n_pnt < 1) {
        std::string msg = "No valid CTA observation has been found in "
                          "observation list, hence no pointing information "
                          "could be extracted. Use the \"usepnt=no\" "
                          "option and specify the pointing explicitly.";
        throw GException::invalid_value(G_GET_MEAN_POINTING, msg);
    }

    // Calculate mean coordinates
    double ra  = xref / double(n_pnt);
    double dec = yref / double(n_pnt);

    // Initialise sky dir
    GSkyDir dir = GSkyDir();
    dir.radec_deg(ra,dec);

    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief Get current resident set size (physical memory use) in Bytes
 *
 * @return Physical memory use in Bytes.
 ***************************************************************************/
size_t ctool::get_current_rss(void)
{
    // Initialize resident set size
    size_t rss = 0;

    // Determine resident set size (architecture dependent)
    // OSX
    #if defined(__APPLE__) && defined(__MACH__)
    #ifdef MACH_TASK_BASIC_INFO
    struct mach_task_basic_info info;
    mach_msg_type_number_t      infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount) == KERN_SUCCESS) {
        rss = (size_t)info.resident_size;
    }
    #else
    struct task_basic_info info;
    mach_msg_type_number_t info_count = TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), TASK_BASIC_INFO,
        (task_info_t)&info, &info_count) == KERN_SUCCESS) {
        rss = (size_t)info.resident_size;
    }
    #endif
    // Linux
    #elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    FILE* fp = NULL;
    if ((fp = fopen( "/proc/self/statm", "r" )) != NULL) {
        if (fscanf( fp, "%*s%ld", &rss ) == 1) {
            rss *= (size_t)sysconf(_SC_PAGESIZE);
        }
        fclose(fp);
    }
    // AIX, BSD, Solaris, and Unknown OS
    #else
    rss = 0;
    #endif

    // Return resident set size
    return rss;
}


/***********************************************************************//**
 * @brief Dumps help text in the console
 ***************************************************************************/
void ctool::provide_help(void) const
{
    // Allocate line buffer
    const int n = 1000; 
    char  line[n];

    // Build help file name
    std::string helpfile = name()+".txt";

    // Get ctools environment variable
    char* ptr = std::getenv("CTOOLS");
    if (ptr == NULL) {
        std::string msg = "CTOOLS environment variable not set, cannot "
                          "display help file. Please set the CTOOLS "
                          "environment variable.";
        throw GException::invalid_value(G_PROVIDE_HELP, msg);
    }

    // If help file exists then display it, otherwise notify that no
    // help is available
    std::string fname = std::string(ptr) + "/share/help/" + helpfile;
    FILE* fptr = fopen(fname.c_str(), "r");
    if (fptr != NULL) {
        while (fgets(line, n, fptr) != NULL) {
            std::cout << std::string(line);
        }
        fclose(fptr);
    }
    else {
        std::cout << "No help available for "+name()+"." << std::endl;
    }

    // Return
    return;
}
