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

    // Initialise run parameters
    m_livetime = 0.0;

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
    m_bkgmdl = m_obs.models();
    m_outmdl.clear();
    m_cube_model = -1;

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
                    log.header3("Observation "+gammalib::str(i));
                }
            }

            // Increment number of observations
            n_observations++;

            // Fill the cube
            fill_cube(obs);

        } // endif: CTA observation found

    } // endfor: looped over observations

    // Re-normalize cube to get units of counts/MeV/s/sr
    if (m_livetime > 0.0) {
    
        // Loop over all bins in background cube and divide the content
        // by the total livetime.
        for (int i = 0; i < m_bkgcube.size(); ++i) {
            GCTAEventBin* bin  = m_bkgcube[i];
            double        rate = bin->counts() / m_livetime;
            bin->counts(rate);
        }

    } // endif: livetime was positive

    // Log results
    if (logTerse()) {
        log << gammalib::parformat("Number of observations");
        log << n_observations << std::endl;
    }

    // Create a background model for output background cube and append
    // that model to the input model in place of the original
    // background model
    GEnergies energies;
    for (int i = 0; i < m_bkgcube.ebins(); ++i) {
        energies.append(m_bkgcube.energy(i));
    }
    GModelSpatialDiffuseCube spatial(m_bkgcube.map(), energies, 1.0);
    GModelSpectralPlaw       spectral(1.0, 0.0, GEnergy(1.0, "TeV"));
    GCTAModelCubeBackground  model(spatial, spectral);
    model.name("ctbkgcube default background model");
    m_cube_model = m_outmdl.size(); // Store the slot number for save()
    m_outmdl.append(model);

    // Log output model
    if (logTerse()) {
        log << std::endl;
        log.header1("Output model");
        log << m_outmdl << std::endl;
    }

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

    // Create energies container from energy boundaries
    GEnergies energies;
    for (int i = 0; i < m_ebounds.size(); ++i) {
        energies.append(m_ebounds.elogmean(i));
    }

    // Get output filenames
    std::string outfile  = (*this)["outcube"].filename();
    std::string outmodel = (*this)["outmodel"].filename();

    // Create empty FITS file
    GFits fits;

    // Write background cube
    m_bkgcube.map().write(fits);

    // Write energies
    energies.write(fits);
    
    // Save FITS file
    fits.saveto(outfile, clobber());

    // Write output models if filename is valid
    if ((outmodel.length() > 0) && (gammalib::tolower(outmodel) != "none")) {

        // Set filename of map cube.
        if (m_cube_model >= 0) {
            GCTAModelCubeBackground* model =
                dynamic_cast<GCTAModelCubeBackground*>(m_outmdl[m_cube_model]);
            if (model != NULL) {
                GModelSpatialDiffuseCube* spatial =
                    dynamic_cast<GModelSpatialDiffuseCube*>(model->spatial());
                if (spatial != NULL) {
                    spatial->filename(outfile);
                }
            }
        }

        // Save output model container
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
    m_bkgcube.clear();
    m_bkgmdl.clear();
    m_outmdl.clear();
    m_ebounds.clear();
    m_livetime   = 0.0;
    m_cube_model = -1;

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
    m_bkgcube    = app.m_bkgcube;
    m_bkgmdl     = app.m_bkgmdl;
    m_outmdl     = app.m_outmdl;
    m_ebounds    = app.m_ebounds;
    m_livetime   = app.m_livetime;
    m_cube_model = app.m_cube_model;

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

    // Get the incube filename
    std::string incube = (*this)["incube"].filename();

    // Check for filename validity
    if ((gammalib::toupper(incube) == "NONE") ||
        (gammalib::strip_whitespace(incube) == "")) {

        // Create an event cube based on task parameters
        m_bkgcube = create_cube(m_obs);

    } // endif: filename was not valid

    // ... otherwise setup the background cube from the counts map
    else {
    
        // Load event cube from filename
        m_bkgcube.load(incube);
        // Loop over all bins in background cube
        for (int i = 0; i < m_bkgcube.size(); ++i) {

            // Get event bin
            GCTAEventBin* bin = m_bkgcube[i];
            bin->counts(0.0);

        }
    
    } // endelse: cube was loaded from file

    // Get energy definition
    m_ebounds = m_bkgcube.ebounds();

    // Read output filenames (if needed)
    if (read_ahead()) {
        m_outcube  = (*this)["outcube"].filename();
        m_outmodel = (*this)["outmodel"].filename();
    }

    // Return
    return;
}

/***********************************************************************//**
 * @brief Generate background cube
 *
 * @param[in] obs Pointer to CTA observation.
 *
 * Fills the background cube with the model value for a CTA observation.
 * This method will not need to load the event file as it gets the event
 * ROI and GTI directly from the observation.
 ***************************************************************************/
void ctbkgcube::fill_cube(GCTAObservation* obs)
{

    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Initialise statistics
        double sum            = 0.0;
        int    n_bins_inside  = 0;
        int    n_bins_outside = 0;

        // Extract region of interest from CTA observation
        GCTARoi roi = obs->roi();

        // Set GTI of actual observations as the GTI of the event cube
        m_bkgcube.gti(obs->gti());

        // Get observation livetime
        double livetime = obs->livetime();

        // Loop over all bins in background cube
        for (int i = 0; i < m_bkgcube.size(); ++i) {

            // Get event bin
            GCTAEventBin* bin = m_bkgcube[i];

            // Continue only if binned in contained in ROI
            if (roi.contains(*bin)) {
            
                // Compute model value for event bin. The model value is
                // given in counts/MeV/s/sr.
                double model = m_bkgmdl.eval(*bin, *obs);

                // Compute number events
                sum += model * bin->size();

                // Multiply by livetime to get the correct weighting for
                // each observation. We divide by the total livetime later
                // to get the background model in units of counts/MeV/s/sr.
                model *= livetime;

                // Add existing number of counts
                model += bin->counts();

                // Store cumulated value (units: counts/MeV/sr)
                bin->counts(model);

                // Increment contained bin counter
                n_bins_inside++;
                
            } // endif: bin was contained in RoI
            
            // ... otherwise count bin as an outsider
            else {
                n_bins_outside++;
            }

        } // endfor: looped over all bins

        // Accumulate livetime
        m_livetime += livetime;

        // Log results
        if (logTerse()) {
            log << gammalib::parformat("Bins within RoI");
            log << n_bins_inside << std::endl;
            log << gammalib::parformat("Bins outside RoI");
            log << n_bins_outside << std::endl;
            log << gammalib::parformat("Background events in cube");
            log << sum << std::endl;
            log << gammalib::parformat("Cube livetime");
            log << livetime << " sec" << std::endl;
        }

    } // endif: observation pointer was not valid

    // Return
    return;
}
