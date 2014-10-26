/***************************************************************************
 *                 ctbutterfly - butterfly calculation tool                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Michael Mayer                                    *
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
 * @file ctbutterfly.cpp
 * @brief Butterfly computation tool
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include <fstream>
#include "ctbutterfly.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                      "ctbutterfly::get_parameters()"
#define G_SAVE                                          "ctbutterfly::save()"

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
ctbutterfly::ctbutterfly(void) : ctool(CTBUTTERFLY_NAME, CTBUTTERFLY_VERSION)
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
ctbutterfly::ctbutterfly(const GObservations& obs) :
             ctool(CTBUTTERFLY_NAME, CTBUTTERFLY_VERSION)
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
ctbutterfly::ctbutterfly(int argc, char *argv[]) :
             ctool(CTBUTTERFLY_NAME, CTBUTTERFLY_VERSION, argc, argv)
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
ctbutterfly::ctbutterfly(const ctbutterfly& app) : ctool(app)
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
ctbutterfly::~ctbutterfly(void)
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
ctbutterfly& ctbutterfly::operator=(const ctbutterfly& app)
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
void ctbutterfly::clear(void)
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
 * @brief Computes the butterfly
 *
 * This method takes the spectral fit and its covariance matrix to compute
 * the confidence band via Gaussian error propagation
 ***************************************************************************/
void ctbutterfly::run(void)
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

    // Get model instance for further computations
    GModels models = m_obs.models();

    // Compute covariance matrix if necessary
    if (m_covariance.size() == 0) {

        // Write header
        if (logTerse()) {
            log << std::endl;
            log.header1("Compute covariance matrix");
        }

        // Evaluate curvature matrix at the actual parameters
        GOptimizerPars            pars       = models.pars();
        GObservations::likelihood likelihood = m_obs.function();
        likelihood.eval(pars);

        // Get covariance matrix by inverting the curvature matrix
        m_covariance = likelihood.curvature()->invert();

    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Generate butterfly");
    }

    // Initialise dummy time to evaluate spectral model
    GTime time = GTime();

    // Initialise vector of gradients
    GVector grad = GVector(m_obs.models().npars());

    // Loop over energy bins positions
    for (int i = 0; i < m_ebounds.size(); ++i) {

        // Initialise number of gradients
        int num_gradient = 0;

        // Get the energy of current bin
        GEnergy energy = m_ebounds.elogmean(i);

        // Header for verbose logging
        if (logExplicit()) {
            std::string msg = "Computing butterfly for bin number "+
                              gammalib::str(i)+" at "+energy.print();
            log << std::endl;
            log.header2(msg);
        }

        // Initialise model flux value
        double model_flux = 0.0;

        // Loop over models
        for (int j = 0; j < models.size(); ++j) {

            // Check wether model is a skymodel
            GModelSky* skymodel = dynamic_cast<GModelSky*>(models[j]);

            // Yes ...
            if (skymodel != NULL) {

                // Skip spatial models
                num_gradient += skymodel->spatial()->size();

                // Get pointer to spectral model
                GModelSpectral* spectral = skymodel->spectral();

                // Set flux value of source of interest
                if (skymodel->name() == m_srcname) {
                    model_flux = spectral->eval_gradients(energy, time);
                } // endif: model was source of interest
                
                // Loop over model parameters, get gradients
                // and assign them to the vector
                for (int k = 0; k < spectral->size(); ++k) {
                    grad[num_gradient] = (*spectral)[k].gradient();
                    num_gradient++;
                } // endfor: looped over spectral parameters

                // Skip temporal models
                num_gradient += skymodel->temporal()->size();

            } // endif: model was sky model

            else {

                // Skip other models (e.g. GModelData instances)
                num_gradient += models[j]->size();

            } // endelse: Model was not GModelSky

        } // endfor: Looped over models

        // Multiply covariance to the gradient vector
        GVector vector = m_covariance * grad;

        // Get the error from the scalar product
        double error = std::sqrt(grad * vector);

        // Store flux, value and energy for saving
        m_fluxes.push_back(model_flux);
        m_energies.push_back(energy.MeV());
        m_errors.push_back(error);

        // Log information
        if (logExplicit()) {
            log << " Flux at ";
            log << energy.print();
            log << ": ";
            log << model_flux;
            log << " +- ";
            log << error;
            log << " ph/cm2/s/MeV";
        }

    } //endfor: loop over energy bin

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save butterfly diagram
 *
 * Saves the butterfly diagram into an ascii file using a column separated
 * value (CSV) format with blanks as separators.
 ***************************************************************************/
void ctbutterfly::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save Butterfly to file");
    }

    // Get output filename
    m_outfile = (*this)["outfile"].filename();

    // Create CSV table with 3 columns
    GCsv table(m_fluxes.size(), 3);

    // Fill CSV table
    for (int i = 0; i < m_fluxes.size(); ++i) {
        table.real(i, 0, m_energies[i]);
        table.real(i, 1, m_fluxes[i]);
        table.real(i, 2, m_errors[i]);
    }

    // Save CSV table
    table.save(m_outfile, " ", clobber());

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
void ctbutterfly::init_members(void)
{
    // Initialise members
    m_infile.clear();
    m_srcname.clear();
    m_outfile.clear();
    m_ebounds.clear();

    // Initialise protected members
    m_obs.clear();
    m_covariance.clear();
    m_energies.clear();
    m_fluxes.clear();
    m_errors.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctbutterfly::copy_members(const ctbutterfly& app)
{
    // Copy attributes
    m_infile  = app.m_infile;
    m_srcname = app.m_srcname;
    m_outfile = app.m_outfile;
    m_ebounds = app.m_ebounds;

    // Copy protected members
    m_obs        = app.m_obs;
    m_covariance = app.m_covariance;
    m_energies   = app.m_energies;
    m_fluxes     = app.m_fluxes;
    m_errors     = app.m_errors;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctbutterfly::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * @exception GException::invalid_value
 *            Test source not found or no RA/DEC parameters found for test
 *            source.
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. Most parameters are only required if no observation exists so
 * far in the observation container. In this case, a single CTA observation
 * will be added to the container, using the definition provided in the
 * parameter file.
 ***************************************************************************/
void ctbutterfly::get_parameters(void)
{
    // If there are no observations in container then add a single CTA
    // observation using the parameters from the parameter file
    if (m_obs.size() == 0) {

        // Allocate CTA observation
        GCTAObservation obs;

        // Get event file name
        std::string filename = (*this)["infile"].filename();

        // Try first to open as FITS file
        try {

            // Load data
            obs.load(filename);

            // Set response
            set_obs_response(&obs);

            // Append observation to container
            m_obs.append(obs);

        }

        // ... otherwise try to open as XML file
        catch (GException::fits_open_error &e) {

            // Load observations from XML file
            m_obs.load(filename);

            // Check if all observations have response information. If
            // not, get the calibration database parameters and set
            // the response properly
            set_response(m_obs);

        } // endcatch: file was an XML file

    } // endif: there was no observation in the container

    // If there is are no models associated with the observations then
    // load now the model definition
    if (m_obs.models().size() == 0) {

        // Get models XML filename
        std::string filename = (*this)["srcmdl"].filename();

        // Setup models for optimizing.
        m_obs.models(GModels(filename));

    } // endif: no models were associated with observations

    // Get name of test source and check container for this name
    m_srcname = (*this)["srcname"].string();
    if (!m_obs.models().contains(m_srcname)) {
        std::string msg = "Source \""+m_srcname+"\" not found in model "
                          "container. Please add a source with that name "
                          "or check for possible typos.";
    	throw GException::invalid_value(G_GET_PARAMETERS, msg);
    }

    // Get energy binning information
    m_ebounds  = get_ebounds();

    // Get matrix file name and load if possible
    std::string matrixfilename = (*this)["matrix"].filename();
    if (matrixfilename != "NONE") {
        std::string msg = "Loading of matrix from file not implemented yet. "
                           "Use filename = \"NONE\" to induce a recomputation of "
                           "the matrix internally.";
        throw GException::feature_not_implemented(G_GET_PARAMETERS, msg);
        // m_covariance.load(matrixfilename);
    }

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outfile = (*this)["outfile"].filename();
    }

    // Return
    return;
}
