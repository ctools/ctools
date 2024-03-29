/***************************************************************************
 *                    cttsmap - TS map calculation tool                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2022 by Michael Mayer                               *
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
 *
 * Constructs an empty cttsmap tool.
 ***************************************************************************/
cttsmap::cttsmap(void) : ctlikelihood(CTTSMAP_NAME, VERSION)
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
 * Constructs cttsmap tool from an observations container.
 ***************************************************************************/
cttsmap::cttsmap(const GObservations& obs) :
         ctlikelihood(CTTSMAP_NAME, VERSION, obs)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Command line constructor
 *
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 *
 * Constructs an instance of the cttsmap tool that will parse user parameters
 * that are provided as command line arguments.
 ***************************************************************************/
cttsmap::cttsmap(int argc, char *argv[]) :
         ctlikelihood(CTTSMAP_NAME, VERSION, argc, argv)
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
cttsmap::cttsmap(const cttsmap& app) : ctlikelihood(app)
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
        this->ctlikelihood::operator=(app);

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
    this->ctlikelihood::free_members();
    this->ctobservation::free_members();
    this->ctool::free_members();

    // Clear base class (needed to conserve tool name and version)
    this->GApplication::clear();

    // Initialise members
    this->ctool::init_members();
    this->ctobservation::init_members();
    this->ctlikelihood::init_members();
    init_members();

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
void cttsmap::process(void)
{
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

    // Set optimizer logger
    if (logExplicit()) {
        static_cast<GOptimizerLM*>(&m_opt)->logger(&log);
    }
    else {
        static_cast<GOptimizerLM*>(&m_opt)->logger(NULL);
    }

    // Store initial models
    GModels models_orig = m_obs.models();

    // Get model instance for further computations
    GModels models = m_obs.models();

    // Remove test source
    models.remove(m_srcname);

    // Fix spatial parameters if requested
    if ((*this)["fix_spat"].boolean()) {

        // Loop over all models
        for (int i = 0; i < models.size(); ++i) {

            // Continue only if model is skymodel
            GModelSky* sky= dynamic_cast<GModelSky*>(models[i]);
            if (sky != NULL) {

                // Fix all spatial parameters
                GModelSpatial* spatial = sky->spatial();
                for (int j = 0; j < spatial->size(); j++) {
                    (*spatial)[j].fix();
                }

            } // endif: there was a sky model

        } // endfor: looped over models

    } // endif: spatial parameter should be fixed

    // Get likelihood of null hypothesis if not given before
    if (m_logL0 == 0.0) {

        // Write header
        log_header1(TERSE, "Compute NULL Hypothesis for TS computation");

        // Compute likelihood without the test source
        m_obs.models(models);
        m_obs.optimize(m_opt);
        m_logL0 = -(m_opt.value());

        // Write input model container into logger
        log_models(VERBOSE, m_obs.models(), "Fitted model");

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
        (*m_testsource)["RA"].value(bincentre.ra_deg());
        (*m_testsource)["DEC"].value(bincentre.dec_deg());
        models.append(*m_testsource);

        // Assign models to observations
        m_obs.models(models);

        // Make sure that no values are cached for test source
        m_obs.remove_response_cache(m_srcname);

        // Optimize observation container
        m_obs.optimize(m_opt);

        // Compute errors if necessary
        if (m_errors) {
            m_obs.errors(m_opt);
        }

        // Get status of optimization
        int status = m_opt.status();

        // Retrieve the Likelihood value
        double logL1 = -(m_opt.value());

        // Compute TS value
        double ts = 2.0 * (logL1 - m_logL0);

        // Log information
        if (logExplicit()) {
            log_value(EXPLICIT, "TS value", ts);
        }
        else if (logNormal()) {
            log_value(NORMAL, "TS value (bin "+gammalib::str(i)+")",
                      gammalib::str(ts)+" ("+bincentre.print()+")");
        }

        // Get test source model instance
        GModels best_fit_model = m_obs.models();
        GModel* testsource     = best_fit_model[m_srcname];

        // Assign values to the maps
        m_tsmap(i) = ts;

        // Extract fitted test source parameters
        for (int j = 0; j < m_mapnames.size(); ++j) {

            // If map name start with "e_" then set the parameter error ...
            if (m_mapnames[j].substr(0,2) == "e_") {

                // Get parameter name by removing the error prefix
                std::string parname = m_mapnames[j].substr(2, m_mapnames[j].size());

                // Get parameter error
                m_maps[j](i) = (*testsource)[parname].error();

            }

            // ... otherwise set the parameter value
            else {
                m_maps[j](i) = (*testsource)[m_mapnames[j]].value();
            }

        } // endfor: looped over all maps

        // Set Fit status of the bin
        m_statusmap(i) = status;

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

    // Save only if filename is valid and TS map is not empty
    if ((*this)["outmap"].is_valid() && !m_tsmap.is_empty()) {

        // Get TS map filename
        m_outmap = (*this)["outmap"].filename();

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

        // Write Fit status map
        m_statusmap.write(fits);
        fits[fits.size()-1]->extname("STATUS MAP");

        // Stamp FITS file
        stamp(fits);

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
    m_errors      = false;

    // Initialise protected members
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
    m_errors      = app.m_errors;

    // Copy protected members
    m_binmin    = app.m_binmin;
    m_binmax    = app.m_binmax;
    m_logL0     = app.m_logL0;
    m_tsmap     = app.m_tsmap;
    m_statusmap = app.m_statusmap;
    m_mapnames  = app.m_mapnames;
    m_maps      = app.m_maps;

    // Clone protected members
    m_testsource = (app.m_testsource != NULL) ? app.m_testsource->clone()
                                              : NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void cttsmap::free_members(void)
{
    // Free test source
    if (m_testsource != NULL) {
        delete m_testsource;
        m_testsource = NULL;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * @exception GException::invalid_value
 *            Test source has no RA/DEC parameters.
 *
 * Get all task parameters from parameter file.
 ***************************************************************************/
void cttsmap::get_parameters(void)
{
    // Setup observations from "inobs" parameter
    setup_observations(m_obs);

    // Set observation statistic
    set_obs_statistic(gammalib::toupper((*this)["statistic"].string()));

    // Setup models from "inmodel" parameter
    setup_models(m_obs, (*this)["srcname"].string());

    // Get name of test source
    m_srcname = (*this)["srcname"].string();

    // Get model and store it as protected member. Note that the clone
    // method will never return a NULL pointer, hence we do not need
    // to check the pointer.
    m_testsource = m_obs.models()[m_srcname]->clone();

    // Check if model has RA and DEC parameters which are necessary to
    // continue
    if (!m_testsource->has_par("RA") || !m_testsource->has_par("DEC")) {
        std::string msg = "Source \""+m_srcname+"\" has no \"RA\" and "
                          "\"DEC\" parameters. Only sources with \"RA\" "
                          "and \"DEC\" parameters can be used as test "
                          "sources.";
        throw GException::invalid_value(G_GET_PARAMETERS, msg);
    }

    // Fix the spatial parameters of the test source
    (*m_testsource)["RA"].fix();
    (*m_testsource)["DEC"].fix();

    // Create sky map based on task parameters
    GSkyMap map = create_map(m_obs);

    // Check if errors should be computed
    m_errors = (*this)["errors"].boolean();

    // Initialise maps from user parameters
    init_maps(map);

    // Read energy dispersion flag
    m_apply_edisp = (*this)["edisp"].boolean();

    // Get optional splitting parameters
    m_binmin  = (*this)["binmin"].integer();
    m_binmax  = (*this)["binmax"].integer();
    m_logL0   = (*this)["logL0"].real();
    m_publish = (*this)["publish"].boolean();

    // Query other parameters
    (*this)["fix_spat"].boolean();

    // If needed later, query output filename now
    if (read_ahead()) {
        (*this)["outmap"].query();
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

    // Create status map and progress map
    m_statusmap = GSkyMap(map);

    // Set status map to invalid value -1 to signal that bin wasn't computed yet
    m_statusmap = -1.0;

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

        // Add parameter errors if requested
        if (m_errors) {

            // Add sky map for fit error of the free parameter
            m_maps.push_back(map);

            // Store name of parameter error
            m_mapnames.push_back("e_" + (*m_testsource)[i].name());

        } // endif: error calculation was requested

    } // endfor: looped over all model parameters

    // Return
    return;
}
