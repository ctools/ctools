/***************************************************************************
 *                  ctmodel - Model cube generation tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2015 by Juergen Knoedlseder                         *
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
 * @file ctmodel.cpp
 * @brief Model cube generation tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctmodel.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                          "ctmodel::get_parameters()"

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
ctmodel::ctmodel(void) : ctool(CTMODEL_NAME, CTMODEL_VERSION)
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
ctmodel::ctmodel(const GObservations& obs) : ctool(CTMODEL_NAME, CTMODEL_VERSION)
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
ctmodel::ctmodel(int argc, char *argv[]) :
         ctool(CTMODEL_NAME, CTMODEL_VERSION, argc, argv)
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
ctmodel::ctmodel(const ctmodel& app) : ctool(app)
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
ctmodel::~ctmodel(void)
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
ctmodel& ctmodel::operator=(const ctmodel& app)
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
void ctmodel::clear(void)
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
 * @brief Generate the model map(s)
 *
 * This method reads the task parameters from the parfile, sets up the
 * observation container, loops over all CTA observations in the container
 * and generates a model map for each CTA observation.
 ***************************************************************************/
void ctmodel::run(void)
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

    // Write models into logger
    if (logTerse()) {
        log << std::endl;
        log.header1("Models");
        log << m_obs.models() << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Generate model cube");
    }

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

            // Fill cube and leave loop if we are binned mode (meaning we only have
            // one binned observation)
            if (m_binned && m_obs.size() == 1) {
                fill_cube(obs);
                break;
            }

            // Skip observation if we don't have an unbinned observation
            if (obs->eventtype() != "EventList") {

                // Log that we skip the this observation
                if (logTerse()) {
                    log << "Warning: Skipping binned observation \""+obs->name()+"\"";
                    log << std::endl;
                }
                continue;
            }


            // Fill the cube
            fill_cube(obs);

            // Leave loop if we are binned mode (meaning we only have
            // one binned observation
            if (m_binned && m_obs.size() == 1) {
                break;
            }

            // Dispose events to free memory if event file exists on disk
            if (obs->eventfile().length() > 0 &&
                gammalib::file_exists(obs->eventfile())) {
                obs->dispose_events();
            }

        } // endif: CTA observation found

    } // endfor: looped over observations

    // Log cube
    if (logTerse()) {
        log << std::endl;
        log.header1("Model cube");
        log << m_cube << std::endl;
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
 * @brief Save model cube
 *
 * Saves the model cube into a FITS file specified using the "outfile"
 * task parameter.
 ***************************************************************************/
void ctmodel::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save cube");
    }

    // Make sure we have the FITS filename
    m_outcube = (*this)["outcube"].filename();

    // Save model cube into FITS file
    m_cube.save(m_outcube, clobber());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set model cube
 *
 * @param[in] cube Model cube.
 *
 * Set model cube and set all cube bins to zero.
 ***************************************************************************/
void ctmodel::cube(const GCTAEventCube& cube)
{
    // Set cube
    m_cube = cube;

    // Set all cube bins to zero
    for (int i = 0; i < m_cube.size(); ++i) {
        m_cube[i]->counts(0.0);
    }

    // Signal that cube has been set
    m_has_cube = true;

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
void ctmodel::init_members(void)
{
    // Initialise members
    m_outcube.clear();
    m_apply_edisp = false;

    // Initialise protected members
    m_obs.clear();
    m_cube.clear();
    m_gti.clear();
    m_has_cube    = false;
    m_append_cube = false;
    m_binned = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctmodel::copy_members(const ctmodel& app)
{
    // Copy attributes
    m_outcube     = app.m_outcube;
    m_apply_edisp = app.m_apply_edisp;

    // Copy protected members
    m_obs         = app.m_obs;
    m_cube        = app.m_cube;
    m_gti         = app.m_gti;
    m_has_cube    = app.m_has_cube;
    m_append_cube = app.m_append_cube;
    m_binned    = app.m_binned;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctmodel::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. The parameters are read in the correct order.
 ***************************************************************************/
void ctmodel::get_parameters(void)
{
    // Reset cube append flag
    m_append_cube = false;

    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {
        m_obs = get_observations();
    }

    // Check if we got excactly one binned CTA observation
    if (m_obs.size() == 1) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

        // Continue only if observation is a CTA observation
        if (obs != NULL) {

            // Check for binned observation
            if (obs->eventtype() == "CountsCube") {

                // Set cube from binned observation
                GCTAEventCube* evtcube = dynamic_cast<GCTAEventCube*>(const_cast<GEvents*>(obs->events()));

                cube(*evtcube);

                // Signal that cube has been set
                m_has_cube = true;

                // Signal that we are in binned mode
                m_binned = true;

            } // endif: observation was binned

        } // endif: observation was CTA

    } // endif: had exactly one observation


    // Read model definition file if required
    if (m_obs.models().size() == 0) {

        // Get model filename
        std::string inmodel = (*this)["inmodel"].filename();

        // Load models from file
        m_obs.models(inmodel);

    } // endif: there were no models

    // Get energy dispersion flag parameters
    m_apply_edisp = (*this)["edisp"].boolean();

    // Optionally get input counts cube for model cube definition
    if (!m_has_cube) {

        // Read cube definition file
        std::string incube = (*this)["incube"].filename();

        // If no cube file has been specified then create a cube from
        // the task parameters ...
        if ((incube == "NONE") ||
            (gammalib::strip_whitespace(incube) == "")) {
            
            // Create cube from scratch
            m_cube = create_cube(m_obs);

        }

        // ... otherwise load the cube from file and reset all bins
        // to zero
        else {

            // Load cube from given file
            m_cube.load(incube);

            // Set all cube bins to zero
            for (int i = 0; i < m_cube.size(); ++i) {
                m_cube[i]->counts(0.0);
            }
        }

        // Signal that cube has been set
        m_has_cube = true;
    }

    // Read optionally output cube filenames
    if (read_ahead()) {
           m_outcube = (*this)["outcube"].filename();
    }

    // If cube should be appended to first observation then do that now.
    // This is a kluge that makes sure that the cube is passed as part
    // of the observation in case that a cube response is used. The kluge
    // is needed because the GCTACubeSourceDiffuse::set method needs to
    // get the full event cube from the observation.
    if (m_append_cube) {

        //TODO: Check that energy boundaries are compatible

        // Attach GTI of observations to model cube
        m_cube.gti(m_obs[0]->events()->gti());
    
        // Attach model cube to observations
        m_obs[0]->events(m_cube);

    } // endif: cube was scheduled for appending

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill model into model cube
 *
 * @param[in] obs CTA observation.
 *
 * Adds the expected number of events for a given observation to the events
 * that are already found in the model cube. The method also updates the
 * GTI of the model cube so that cube GTI is a list of the GTIs of all
 * observations that were used to generate the model cube.
 ***************************************************************************/
void ctmodel::fill_cube(const GCTAObservation* obs)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Get energy boundaries and GTI references for observation
        const GEbounds& ebounds = obs->events()->ebounds();
        const GGti&     gti     = obs->events()->gti();

        // Initialise statistics
        double sum              = 0.0;
        int    num_outside_ebds = 0;

        // Setup cube GTIs for this observation
        m_cube.gti(obs->events()->gti());

        // Loop over all cube bins
        for (int i = 0; i < m_cube.size(); ++i) {

            // Get cube bin
            GCTAEventBin* bin = m_cube[i];

            // Skip bin if it is outside the energy range of the observation
            if (!ebounds.contains(bin->energy())) {
                num_outside_ebds++;
                continue;
            }

            // Get actual bin value
            double value = bin->counts();
            
            // Compute model value for cube bin
            double model = m_obs.models().eval(*bin, *obs) * bin->size();
            //double model = m_obs.models().eval(*bin, *ptr) * bin->size();

            // Add model to actual value
            value += model;
            sum   += model;

            // Store value
            bin->counts(value);

        } // endfor: looped over all cube bins

        // Append GTIs of observation to list of GTIs
        m_gti.extend(gti);

        // Update GTIs
        m_cube.gti(m_gti);

        // Log results
        if (logTerse()) {
            log << gammalib::parformat("Model events in cube");
            log << sum << std::endl;
            log << gammalib::parformat("Bins outside energy range");
            log << num_outside_ebds << std::endl;
        }

        // Log cube
        if (logExplicit()) {
            log.header2("Model cube");
            log << m_cube << std::endl;
        }

    } // endif: observation was valid

    // Return
    return;
}

