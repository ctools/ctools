/***************************************************************************
 *                 ctexpcube - Exposure cube generation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file ctexpcube.cpp
 * @brief Exposure cube generation tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctexpcube.hpp"
#include "GTools.hpp"
#include "GWcs.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                        "ctexpcube::get_parameters()"
#define G_SET_FROM_CNTMAP          "ctexpcube::set_from_cntmap(std::string&)"

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
ctexpcube::ctexpcube(void) : ctool(CTEXPCUBE_NAME, CTEXPCUBE_VERSION)
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
ctexpcube::ctexpcube(const GObservations& obs) :
           ctool(CTEXPCUBE_NAME, CTEXPCUBE_VERSION)
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
ctexpcube::ctexpcube(int argc, char *argv[]) :
           ctool(CTEXPCUBE_NAME, CTEXPCUBE_VERSION, argc, argv)
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
 * @param[in] app ctexpcube application.
 ***************************************************************************/
ctexpcube::ctexpcube(const ctexpcube& app) : ctool(app)
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
ctexpcube::~ctexpcube(void)
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
 * @param[in] app ctexpcube application.
 * @return Returns ctexpcube application.
 ***************************************************************************/
ctexpcube& ctexpcube::operator=(const ctexpcube& app)
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
 * @brief Clear instance
 ***************************************************************************/
void ctexpcube::clear(void)
{
    // Free members
    free_members();
    this->ctool::free_members();
    this->GApplication::free_members();

    // Initialise members
    this->GApplication::init_members();
    this->ctool::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Execute application
 *
 * This is the main execution method of the ctexpcube class. It is invoked
 * when the executable is called from command line. The method generates
 * the exposure cube and saves the result.
 ***************************************************************************/
void ctexpcube::execute(void)
{
    // Signal that some parameters should be read ahead
    m_read_ahead = true;

    // Create the exposure cube
    run();

    // Save the exposure cube into the output FITS file
    save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generate the exposure cube(s).
 *
 * This method reads the task parameters from the parfile, sets up the
 * observation container, loops over all CTA observations in the container
 * and generates an exposure cube from the CTA observations.
 ***************************************************************************/
void ctexpcube::run(void)
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
        log.header1("Generate exposure cube");
    }

    // Fill exposure
    m_expcube.fill(m_obs);

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
 * @brief Save exposure cube
 ***************************************************************************/
void ctexpcube::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save exposure cube");
    }

    // Get output filename
    m_outfile = (*this)["outfile"].filename();

    // Save exposure cube
    m_expcube.save(m_outfile, clobber());

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
void ctexpcube::init_members(void)
{
    // Initialise members
    m_outfile.clear();
    m_apply_edisp = false;

    // Initialise protected members
    m_read_ahead = false;
    m_obs.clear();
    m_expcube.clear();
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
void ctexpcube::copy_members(const ctexpcube& app)
{
    // Copy attributes
    m_outfile     = app.m_outfile;
    m_apply_edisp = app.m_apply_edisp;

    // Copy protected members
    m_read_ahead = app.m_read_ahead;
    m_obs        = app.m_obs;
    m_expcube    = app.m_expcube;
    m_ebounds    = app.m_ebounds;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctexpcube::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. The parameters are read in the correct order.
 *
 * @todo Setup exposure cube from counts map
 ***************************************************************************/
void ctexpcube::get_parameters(void)
{
    // If we do not have any observations in the container then get an
    // input file name or observation descriptor file
    if (m_obs.size() == 0) {
        get_obs();
    }

    // For all observations that have no response, set the response
    // from the task parameters
    set_response(m_obs);

    // Read energy dispersion flag
    m_apply_edisp = (*this)["edisp"].boolean();

    // If no counts map is specified then setup the exposure cube from
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
        m_ebounds = get_ebounds();

        // Define exposure cube
        m_expcube = GCTAExposure(wcs, coordsys, xref, yref,
                                 -binsz, binsz, nxpix, nypix,
                                 m_ebounds);
    }

    // ... otherwise setup the exposure cube from the counts map
    else {
    
        // Set exposure cube from counts map
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
void ctexpcube::get_obs(void)
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
 * @brief Set exposure cube definition from counts map
 *
 * @exception GException::invalid_value
 *            Invalid counts map projection or invalid events encountered.
 *
 * Set exposure cube definition from counts map.
 ***************************************************************************/
void ctexpcube::set_from_cntmap(const std::string& filename)
{
    // Allocate CTA observation
    GCTAObservation obs;
        
    // Load counts map in CTA observation
    obs.load(filename);

    // Set exposure cube from counts map
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

            // Define exposure cube
            m_expcube = GCTAExposure(proj, coordsys, xref, yref,
                                     dx, dy, nx, ny,
                                     m_ebounds);
        
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
