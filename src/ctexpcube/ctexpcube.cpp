/***************************************************************************
 *                   ctexpcube - CTA exposure cube tool                    *
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
 * @brief CTA exposure cube tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctexpcube.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_EBOUNDS                              "ctexpcube::get_ebounds()"


#define G_SETUP_OBS                                    "ctexpcube::setup_obs()"
#define G_MODEL_MAP                    "ctexpcube::model_map(GCTAObservation*)"

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
ctexpcube::ctexpcube(void) : GApplication(CTEXPCUBE_NAME, CTEXPCUBE_VERSION)
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
 * This method creates an instance of the class by copying an existing
 * observations container.
 ***************************************************************************/
ctexpcube::ctexpcube(GObservations obs) : GApplication(CTEXPCUBE_NAME, CTEXPCUBE_VERSION)
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
           GApplication(CTEXPCUBE_NAME, CTEXPCUBE_VERSION, argc, argv)
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
ctexpcube::ctexpcube(const ctexpcube& app) : GApplication(app)
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
void ctexpcube::clear(void)
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
 * @brief Generate the model map(s)
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
            save_edisp[i] = obs->response().apply_edisp();
            obs->response().apply_edisp(m_apply_edisp);
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

    // Initialise observation counter
    int n_observations = 0;

    // Loop over all observations in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Continue only if observation is a CTA observation
        if (obs != NULL) {

            // Write header for observation
            if (logTerse()) {
                if (obs->name().length() > 1) {
                    log.header3("Observation "+obs->name());
                }
                else {
                    log.header3("Observation");
                }
            }

            // Increment number of observations
            n_observations++;

            // Generate exposure cube
            //model_map(obs, m_obs.models());

        } // endif: CTA observation found

    } // endfor: looped over observations

    // Restore energy dispersion flag for all CTA observations
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            obs->response().apply_edisp(save_edisp[i]);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save exposure cube
 *
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


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. The parameters are read in the correct order.
 ***************************************************************************/
void ctexpcube::get_parameters(void)
{
    // If we do not have any observations in the container then get an
    // input file name or observation descriptor file
    if (m_obs.size() == 0) {
        m_infile = (*this)["infile"].filename();
    }

    // Read response parameters
    m_caldb       = (*this)["caldb"].string();
    m_irf         = (*this)["irf"].string();
    m_apply_edisp = (*this)["edisp"].boolean();

    // If no counts map is specified then setup the exposure cube from
    // the user parameters
    m_cntmap = (*this)["cntmap"].filename();
    if ((m_cntmap == "NONE") || (gammalib::strip_whitespace(m_cntmap) == "")) {
    
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

        // Define exposure cube
        m_expcube = GCTAExposure(wcs, coordsys, xref, yref,
                                 -binsz, binsz, nxpix, nypix,
                                 m_ebounds);
    }

    // ... otherwise setup the exposure cube from the counts map
    else {
    }

    // Read output filename (if needed)
    if (m_read_ahead) {
        m_outfile = (*this)["outfile"].filename();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get energy boundaries from parameters
 *
 * Get the energy boundaries from the user parameters.
 ***************************************************************************/
void ctexpcube::get_ebounds(void)
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
    m_infile.clear();
    m_outfile.clear();
    m_caldb.clear();
    m_irf.clear();
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
    m_infile      = app.m_infile;
    m_outfile     = app.m_outfile;
    m_caldb       = app.m_caldb;
    m_irf         = app.m_irf;
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
    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}
