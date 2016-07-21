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
#define G_GET_PARAMETERS                        "ctbkgcube::get_parameters()"

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
ctbkgcube::ctbkgcube(void) : ctool(CTBKGCUBE_NAME, CTBKGCUBE_VERSION)
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
           ctool(CTBKGCUBE_NAME, CTBKGCUBE_VERSION)
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
ctbkgcube::ctbkgcube(int argc, char *argv[]) :
           ctool(CTBKGCUBE_NAME, CTBKGCUBE_VERSION, argc, argv)
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
ctbkgcube::ctbkgcube(const ctbkgcube& app) : ctool(app)
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
void ctbkgcube::clear(void)
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

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Warn if there are not enough energy bins
    int    n          = m_background.energies().size();
    double logEmin    = std::log10(m_background.energies()[0].TeV());
    double logEmax    = std::log10(m_background.energies()[n-1].TeV());
    int    nrequired  = int((logEmax - logEmin) * 25.0);
    if ((n < nrequired) && logTerse()) {
        log << std::endl;
        log << "WARNING: ";
        log << "Only " << n-1 << " energy bins have been requested. This may be "
               "too few energy bins.";
        log << std::endl << "         ";
        log << "At least 25 bins per decade in energy are recommended, which "
               "in your case";
        log << std::endl << "         ";
        log << "would be " << nrequired-1 << " bins. Consider increasing the "
               "number of energy bins.";
        log << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Initialise background cube");
    }

    // Initialise exposure cube
    init_cube();

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        log.header1(gammalib::number("Observation",m_obs.size()));
        log << m_obs.print(m_chatter) << std::endl;
    }

    // Log input model
    if (logTerse()) {
        log << std::endl;
        log.header1("Input model");
        log << m_obs.models().print(m_chatter) << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Prepare model");
    }

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

        // Do we have a CTA, HESS, MAGIC or VERITAS specific model?
        if (m_bkgmdl[i]->is_valid("CTA", "") ||
            m_bkgmdl[i]->is_valid("HESS", "") ||
            m_bkgmdl[i]->is_valid("MAGIC", "") ||
            m_bkgmdl[i]->is_valid("VERITAS", "")) {

            // Do we have a background model?
            if (dynamic_cast<GModelData*>(m_bkgmdl[i]) != NULL) {

                // ... then keep model
                remove = false;

            } // endif: had a background model
        } // endif: had a CTA, HESS, MAGIC or VERITAS model

        // Log results
        if (logTerse()) {
            if (remove) {
                log << gammalib::parformat("Remove model");
            }
            else {
                log << gammalib::parformat("Keep model");
            }
            log << m_bkgmdl[i]->name();
            log << " ";
            log << m_bkgmdl[i]->type();
            log << "(";
            log << m_bkgmdl[i]->instruments();
            log << ")";
            log << std::endl;
        }

        // If removal is requested, append model to output container and
        // remove it from the background model container ...
        if (remove) {
            m_outmdl.append(*(m_bkgmdl[i]));
            m_bkgmdl.remove(i);
        }

    } // endfor: looped over all background models

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Generate background cube");
    }

    // Assign background models to container
    m_obs.models(m_bkgmdl);

    // Fill background cube from observations
    m_background.fill(m_obs, &log);

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

    // Log background cube
    if (logTerse()) {
        log << m_background.print(m_chatter) << std::endl;
    }

    // Log output model
    if (logTerse()) {
        log << std::endl;
        log.header1("Output model");
        log << m_outmdl << std::endl;
    }

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
    if (logTerse()) {
        log << std::endl;
        log.header1("Save background cube");
    }

    // Get background cube and model file names
    m_outcube  = (*this)["outcube"].filename();
    m_outmodel = (*this)["outmodel"].filename();

    // Save only if filename is non-empty
    if (!m_outcube.is_empty()) {

        // Log filename
        if (logTerse()) {
            log << "Save background cube into file \""+m_outcube+"\".";
            log << std::endl;
        }

        // Save background cube
        m_background.save(m_outcube, clobber());
    
    }

    // Write output models if filename is valid
    if ((!m_outmodel.is_empty()) &&
        (gammalib::tolower(m_outmodel.url()) != "none")) {

        // Log filename
        if (logTerse()) {
            log << "Save model into file \""+m_outmodel.url()+"\".";
            log << std::endl;
        }

        // Save output model for stacked analyses
        m_outmdl.save(m_outmodel.url());

    }

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
    if (logTerse()) {
        log << std::endl;
        log.header1("Publish background cube");
    }

    // Set default name is user name is empty
    std::string user_name(name);
    if (user_name.empty()) {
        user_name = CTBKGCUBE_NAME;
    }

    // Log filename
    if (logTerse()) {
        log << "Publish \""+user_name+"\" background cube." << std::endl;
    }

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
    m_obs.clear();
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
    m_obs        = app.m_obs;
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
 * @exception GException::invalid_value
 *            No background model definition XML file specified.
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. The parameters are read in the correct order.
 ***************************************************************************/
void ctbkgcube::get_parameters(void)
{
    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {

        // Throw exception if no input observation file is given
        require_inobs(G_GET_PARAMETERS);

        // Throw exception if counts cube is given
        require_inobs_nocube(G_GET_PARAMETERS);

        // Build observation container
        m_obs = get_observations();

    } // endif: there was no observation in the container

    // Get the incube filename
    std::string incube = (*this)["incube"].filename();

    // Check for filename validity
    if ((gammalib::toupper(incube) == "NONE") ||
        (gammalib::strip_whitespace(incube) == "")) {

        // Create an event cube based on task parameters
        GCTAEventCube cube = create_cube(m_obs);

        // Define background cube
        m_background = GCTACubeBackground(cube);

    } // endif: filename was not valid

    // ... otherwise setup the background cube from the counts map
    else {

        // Load event cube from filename
        GCTAEventCube cube(incube);

        // Define background cube from file
        m_background = GCTACubeBackground(cube);
    
    } // endelse: cube was loaded from file

    // If there are no models associated with the observations then load now
    // the model definition from the XML file
    if (m_obs.models().size() == 0) {
        if ((*this)["inmodel"].is_undefined()) {
            std::string msg = "No model definition XML file specified. "
                              "Please set the \"inmodel\" parameter to the "
                              "XML file that contains the background model "
                              "definition.";
            throw GException::invalid_value(G_GET_PARAMETERS, msg);
        }
        std::string inmodel = (*this)["inmodel"].filename();
        GModels     models(inmodel);
        m_obs.models(models);
    }

    // Get remaining parameters
    m_addbounds = (*this)["addbounds"].boolean();
    m_publish   = (*this)["publish"].boolean();
    m_chatter   = static_cast<GChatter>((*this)["chatter"].integer());

    // Read output filenames (if needed)
    if (read_ahead()) {
        m_outcube  = (*this)["outcube"].filename();
        m_outmodel = (*this)["outmodel"].filename();
    }

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

        // Set logger
        GLog* logger = NULL;
        if (logTerse()) {
            logger = &log;
        }
    
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
            energies = insert_energy_boundaries(energies, *cta, logger);

        } // endfor: looped over all observations

       } // endif: energy bin insertion requested

    // Setup background cube
    m_background = GCTACubeBackground(wcs, coords, x, y, dx, dy, nx, ny, energies);

    // Log background cube
    if (logTerse()) {
        log << m_background.print(m_chatter) << std::endl;
    }

    // Return
    return;
};
