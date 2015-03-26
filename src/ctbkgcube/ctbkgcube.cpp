/***************************************************************************
 *               ctbkgcube - Background cube generation tool               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Chia-Chun Lu                                *
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
#define G_FILL_CUBE                  "ctbkgcube::fill_cube(GCTAObservation*)"

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

    // Log input model
    if (logTerse()) {
        log << std::endl;
        log.header1("Input model");
        log << m_obs.models() << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Prepare model");
    }

    // Copy models from observation container and reset output model
    // container
    GModels models_orig = m_obs.models();
    m_bkgmdl = m_obs.models();
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
    m_background.fill(m_obs);

    // Create a background model for the output background cube and append
    // that model to the input model in place of the original
    // background models
    // Todo: We might think of creating the spectral model via user parameter
    GCTAModelCubeBackground model(GModelSpectralConst(1.0));

    // Set model name
    model.name("BackgroundModel");

    // Set model instrument
    // Todo: Account for possibility to have observations from different
    // IACTs in the same container
    model.instruments("CTA,HESS,MAGIC,VERITAS");

    // Set model identifier to "0". This is coupled to ctbin, where the binned output
    // observation id is set to 0.
    model.ids("0");

    // Append model to output container
    m_outmdl.append(model);

    // Log output model
    if (logTerse()) {
        log << std::endl;
        log.header1("Output model");
        log << m_outmdl << std::endl;
    }

    // Recover original models
    m_obs.models(models_orig);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save background cube
 *
 * Save the background cube into the file specified by the outfile parameter.
 ***************************************************************************/
void ctbkgcube::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save background cube");
    }


    // Get output filenames
    std::string outfile  = (*this)["outcube"].filename();
    std::string outmodel = (*this)["outmodel"].filename();

    // Save background cube
    m_background.save(outfile, clobber());

    // Write output models if filename is valid
    if ((outmodel.length() > 0) && (gammalib::tolower(outmodel) != "none")) {

        // Save output model for binned analyses
        m_outmdl.save(outmodel);
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
    m_outcube.clear();
    m_outmodel.clear();
    m_obs.clear();
    m_background.clear();
    m_bkgmdl.clear();
    m_outmdl.clear();
    m_ebounds.clear();

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
    // Copy members
    m_outmodel   = app.m_outmodel;
    m_outcube    = app.m_outcube;
    m_obs        = app.m_obs;
    m_background = app.m_background;
    m_bkgmdl     = app.m_bkgmdl;
    m_outmdl     = app.m_outmdl;
    m_ebounds    = app.m_ebounds;

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

        // Define background cube from file
        m_background = GCTACubeBackground(incube);
    
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

//    // Get energy definition
//    m_ebounds = m_background.ebounds();

    // Read output filenames (if needed)
    if (read_ahead()) {
        m_outcube  = (*this)["outcube"].filename();
        m_outmodel = (*this)["outmodel"].filename();
    }

    // Return
    return;
}
