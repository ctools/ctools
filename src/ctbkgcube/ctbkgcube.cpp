/***************************************************************************
 *               ctbkgcube - Background cube generation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Chia-Chun Lu                                *
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
 * @brief Background cube generation tool definition
 * @author Chia-Chun Lu
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctbkgcube.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
ctbkgcube::ctbkgcube(void) : ctobservation(CTBKGCUBE_NAME, CTBKGCUBE_VERSION)
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
ctbkgcube::ctbkgcube(const GObservations& obs) :
		ctobservation(CTBKGCUBE_NAME, CTBKGCUBE_VERSION, obs)
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
 ***************************************************************************/
ctbkgcube::ctbkgcube(int argc, char *argv[]) :
		ctobservation(CTBKGCUBE_NAME, CTBKGCUBE_VERSION, argc, argv)
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
ctbkgcube::ctbkgcube(const ctbkgcube& app) : ctobservation(app)
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
 * @param[in] app Application.
 * @return Returns application.
 ***************************************************************************/
ctbkgcube& ctbkgcube::operator=(const ctbkgcube& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->ctobservation::operator=(app);

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
 * @brief Clear ctbkgcube tool
 *
 * Clears ctbkgcube tool.
 ***************************************************************************/
void ctbkgcube::clear(void)
{
    // Free members
    free_members();
    this->ctobservation::free_members();
    this->ctool::free_members();

    // Clear base class (needed to conserve tool name and version)
    this->GApplication::clear();

    // Initialise members
    this->ctool::init_members();
    this->ctobservation::init_members();
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generate the background cube(s).
 *
 * This method reads the task parameters from the parfile, sets up the
 * observation container, loops over all CTA observations in the container
 * and generates a background cube from the CTA observations.
 ***************************************************************************/
void ctbkgcube::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Warn if there are not enough energy bins
    log_string(TERSE, warn_too_few_energies(m_background.energies()));

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write input model container into logger
    log_models(NORMAL, m_obs.models(), "Input model");

    // Initialise exposure cube
    init_cube();

    // Write header
    log_header1(TERSE, "Prepare model");

    // Copy models from observation container and reset output model
    // container
    GModels models_orig = m_obs.models();
    m_bkgmdl            = m_obs.models();
    m_outmdl.clear();

    // Remove all models that are not CTA background models from the
    // container and put all removed components in the output
    // container
    int num = m_bkgmdl.size();
    for (int i = num-1; i >= 0; --i) {

        // Flag removal
        bool remove = true;

        // If we have a data space model with "GCTA" classname then we
        // have a CTA background model and we want to keep the model
        if ((dynamic_cast<GModelData*>(m_bkgmdl[i]) != NULL) &&
            (m_bkgmdl[i]->classname().substr(0,4) == "GCTA")) {
            remove = false;
        }

        // Log model removal or keeping
        std::string what  = (remove) ? "Remove model" : "Keep model";
        std::string value = m_bkgmdl[i]->name()+" "+m_bkgmdl[i]->type()+"("+
                            m_bkgmdl[i]->instruments()+")";
        log_value(NORMAL, what, value);

        // If removal is requested, append model to output container and
        // remove it from the background model container ...
        if (remove) {
            m_outmdl.append(*(m_bkgmdl[i]));
            m_bkgmdl.remove(i);
        }

    } // endfor: looped over all background models

    // Write header
    log_header1(TERSE, "Generate background cube");

    // Assign background models to container
    m_obs.models(m_bkgmdl);

    // Set pointer to logger dependent on chattiness
    GLog* logger = (logNormal()) ? &log : NULL;

    // Fill background cube from observations
    m_background.fill(m_obs, logger);

    // Create a background model for the output background cube and append
    // that model to the input model in place of the original
    // background models
    // TODO: We might think of creating the spectral model via user parameter
    GModelSpectralPlaw spectral(1.0, 0.0, GEnergy(1.0, "TeV"));
    spectral["Prefactor"].range(0.01, 100.0);
    spectral["Index"].range(-5.0, 5.0);
    GCTAModelCubeBackground model(spectral);

    // Set model name
    model.name("BackgroundModel");

    // Set model instrument
    // TODO: Account for possibility to have observations from different
    // IACTs in the same container
    model.instruments("CTA,HESS,MAGIC,VERITAS");

    // Append model to output container
    m_outmdl.append(model);

    // Write background cube into logger
    log_string(NORMAL, m_background.print(m_chatter));

    // Write input model container into logger
    log_models(NORMAL, m_outmdl, "Output model");

    // Recover original models
    m_obs.models(models_orig);

    // Optionally publish background cube
    if (m_publish) {
        publish();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save background cube
 *
 * Save the background cube into a FITS file and the output model into an
 * XML file.
 ***************************************************************************/
void ctbkgcube::save(void)
{
    // Write header
    log_header1(TERSE, "Save background cube");

    // Get background cube and model file names
    m_outcube  = (*this)["outcube"].filename();
    m_outmodel = (*this)["outmodel"].filename();

    // Determine whether a model definition file should be written
    bool write_moddef = is_valid_filename(m_outmodel);

    // Save background cube if filename and the background cube are not empty
    if (!m_outcube.is_empty() && !m_background.cube().is_empty()) {
        m_background.save(m_outcube, clobber());
    }

    // Save model definition file if requested
    if (write_moddef) {
        m_outmdl.save(m_outmodel.url());
    }

    // Write into logger what has been done
    std::string fname = (m_outcube.is_empty()) ? "NONE" : m_outcube.url();
    std::string mname = (write_moddef) ? m_outmodel.url() : "NONE";
    if (m_background.cube().is_empty()) {
        fname.append(" (cube is empty, no file created)");
    }
    log_value(NORMAL, "Background cube file", fname);
    log_value(NORMAL, "Model definition file", mname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Publish background cube
 *
 * @param[in] name Background cube name.
 ***************************************************************************/
void ctbkgcube::publish(const std::string& name)
{
    // Write header
    log_header1(TERSE, "Publish background cube");

    // Set default name if user name is empty
    std::string user_name(name);
    if (user_name.empty()) {
        user_name = CTBKGCUBE_NAME;
    }

    // Log background cube name
    log_value(NORMAL, "Publish background cube", user_name);

    // Publish exposure cube
    m_background.cube().publish(user_name);

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
    // Initialise user parameters
    m_outcube.clear();
    m_outmodel.clear();
    m_addbounds = true;
    m_publish   = false;
    m_chatter   = static_cast<GChatter>(2);

    // Initialise protected members
    m_background.clear();
    m_bkgmdl.clear();
    m_outmdl.clear();

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
    // Copy user parameters
    m_outmodel  = app.m_outmodel;
    m_outcube   = app.m_outcube;
    m_addbounds = app.m_addbounds;
    m_publish   = app.m_publish;
    m_chatter   = app.m_chatter;

    // Copy protected members
    m_background = app.m_background;
    m_bkgmdl     = app.m_bkgmdl;
    m_outmdl     = app.m_outmdl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctbkgcube::free_members(void)
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
void ctbkgcube::get_parameters(void)
{
    // Setup observations from "inobs" parameter. Require response and accept
    // event lists, but do not accept counts cubes.
    setup_observations(m_obs, true, true, false);

    // Get the incube filename
    std::string incube = (*this)["incube"].filename();

    // If the "incube" file name is valid then setup the background cube from
    // the counts cube. Otherwise create a counts cube from the user
    // parameters
    GCTAEventCube cube = is_valid_filename(incube) ? GCTAEventCube(incube)
                                                   : create_cube(m_obs);
    
    // Define background cube
    m_background = GCTACubeBackground(cube);

    // Setup models from "inmodel" parameter
    setup_models(m_obs);

    // Get remaining parameters
    m_addbounds = (*this)["addbounds"].boolean();
    m_publish   = (*this)["publish"].boolean();
    m_chatter   = static_cast<GChatter>((*this)["chatter"].integer());

    // Read output filenames (if needed)
    if (read_ahead()) {
        m_outcube  = (*this)["outcube"].filename();
        m_outmodel = (*this)["outmodel"].filename();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise background cube
 *
 * Initialise the background cube.
 ***************************************************************************/
void ctbkgcube::init_cube(void)
{
    // Write header
    log_header1(TERSE, "Initialise background cube");

    // Extract exposure cube definition
    const GWcs* proj   = static_cast<const GWcs*>(m_background.cube().projection());
    std::string wcs    = m_background.cube().projection()->code();
    std::string coords = m_background.cube().projection()->coordsys();
    double      x      = proj->crval(0);
    double      y      = proj->crval(1);
    double      dx     = proj->cdelt(0);
    double      dy     = proj->cdelt(1);
    int         nx     = m_background.cube().nx();
    int         ny     = m_background.cube().ny();

    // Extract energies
    GEnergies energies = m_background.energies();

    // If requested, insert energies at all event list energy boundaries
    if (m_addbounds) {

        // Loop over all observations
        for (int i = 0; i < m_obs.size(); ++i) {
    
            // Get observation and continue only if it is a CTA observation
            const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>
                                         (m_obs[i]);

            // Skip observation if it's not a CTA observation
            if (cta == NULL) {
                continue;
            }

            // Skip observation if it does not contain an event list
            if (cta->eventtype() != "EventList") {
                continue;
            }

            // Insert energy boundaries
            energies = insert_energy_boundaries(energies, *cta);

        } // endfor: looped over all observations

       } // endif: energy bin insertion requested

    // Setup background cube
    m_background = GCTACubeBackground(wcs, coords, x, y, dx, dy, nx, ny, energies);

    // Write background cube into logger
    log_string(NORMAL, m_background.print(m_chatter));

    // Return
    return;
}
