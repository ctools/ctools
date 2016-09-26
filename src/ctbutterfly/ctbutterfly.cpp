/***************************************************************************
 *                 ctbutterfly - butterfly calculation tool                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Michael Mayer                               *
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
#define G_CHECK_MODEL                            "ctbutterfly::check_model()"

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
ctbutterfly::ctbutterfly(void) : ctlikelihood(CTBUTTERFLY_NAME, CTBUTTERFLY_VERSION)
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
             ctlikelihood(CTBUTTERFLY_NAME, CTBUTTERFLY_VERSION, obs)
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
ctbutterfly::ctbutterfly(int argc, char *argv[]) :
             ctlikelihood(CTBUTTERFLY_NAME, CTBUTTERFLY_VERSION, argc, argv)
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
ctbutterfly::ctbutterfly(const ctbutterfly& app) : ctlikelihood(app)
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
 * @brief Clear ctbutterfly tool
 *
 * Clears ctbutterfly tool.
 ***************************************************************************/
void ctbutterfly::clear(void)
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

    // Write header into logger
    log_header();

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

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // If fit is selected then do maximum likelihood fit
    if (m_fit) {

        // Write header into logger
        log_header1(TERSE, "Compute best-fit likelihood");

        // Optimize and compute errors
        m_obs.optimize(m_opt);
        m_obs.errors(m_opt);

        // Get covariance matrix by inverting the curvature matrix. Note that
        // the curvature matrix is for the unscaled parameters, hence we have
        // to scale them before we can proceed.
        GObservations::likelihood likelihood = m_obs.function();
        m_covariance = likelihood.curvature()->invert();

        // Write optimizer, optimised model and covariance matrix into logger
        log_string(NORMAL, m_opt.print(m_chatter));
        log_string(NORMAL, m_obs.models().print(m_chatter));
        log_string(EXPLICIT, m_covariance.print(m_chatter));

    } // endif: fit selected

    // Get model instance for further computations
    GModels models = m_obs.models();

    // If no fit is selected and no covariance matrix has been specified
    // then compute now the covariance matrix
    if (!m_fit && m_covariance.size() == 0) {

        // Write header into logger
        log_header1(TERSE, "Compute covariance matrix");

        // Evaluate curvature matrix at the actual parameters
        GOptimizerPars            pars       = models.pars();
        GObservations::likelihood likelihood = m_obs.function();
        likelihood.eval(pars);

        // Get covariance matrix by inverting the curvature matrix. Note that
        // the curvature matrix is for the unscaled parameters, hence we have
        // to scale them before we can proceed.
        m_covariance = likelihood.curvature()->invert();

        // Write covariance matrix into logger
        log_string(EXPLICIT, m_covariance.print(m_chatter));

    } // endif: covariance computation needed

    // Write header into logger
    log_header1(TERSE, "Prepare butterfly computation");

    // Find parameter indices. We do this by the memory location of the
    // model parameters.
    GModelSky* skymodel = dynamic_cast<GModelSky*>(models[m_srcname]);
    GModelPar* par_prefactor = &(*skymodel)["Prefactor"];
    GModelPar* par_index     = &(*skymodel)["Index"];
    int        inx_prefactor = -1;
    int        inx_index     = -1;
    for (int i = 0; i < models.npars(); ++i) {
        if (models.pars()[i] == par_prefactor) {
            inx_prefactor = i;
        }
        if (models.pars()[i] == par_index) {
            inx_index = i;
        }
    }

    // Store model values and scales
    double  prefactor_mean  = (*skymodel)["Prefactor"].value();
    double  index_mean      = (*skymodel)["Index"].value();
    double  pivot_mean      = (*skymodel)["PivotEnergy"].value();
    double  prefactor_scale = (*skymodel)["Prefactor"].scale();
    double  index_scale     = (*skymodel)["Index"].scale();
    GEnergy pivot(pivot_mean, "MeV");

    // Write parameter indices into logger
    log_value(NORMAL, "Prefactor", gammalib::str(prefactor_mean)+
                      " ph/cm2/s/MeV ("+gammalib::str(inx_prefactor)+")");
    log_value(NORMAL, "Index", gammalib::str(index_mean)+
                      " ("+gammalib::str(inx_index)+")");
    log_value(NORMAL, "Pivot energy", pivot.print());

    // Extract covariance elements
    double cov_pp = m_covariance(inx_prefactor,inx_prefactor);
    double cov_pg = m_covariance(inx_prefactor,inx_index);
    double cov_gg = m_covariance(inx_index,inx_index);

    // Compute eigenvectors and eigenvalues
    double  lambda1;
    double  lambda2;
    GVector vector1(2);
    GVector vector2(2);
    eigenvectors(cov_pp, cov_pg, cov_pg, cov_gg,
                 &lambda1, &lambda2, &vector1, &vector2);

    // Write eigenvectors and eigenvalues into logger
    log_value(NORMAL, "Eigenvalue 1", lambda1);
    log_value(NORMAL, "Eigenvalue 2", lambda2);
    log_value(NORMAL, "Eigenvector 1", vector1.print());
    log_value(NORMAL, "Eigenvector 2", vector2.print());

    // Confidence scaling
    double sigma = gammalib::erfinv(m_confidence) * gammalib::sqrt_two;
    double scale = (sigma*sigma);

    // Write confidence level information into logger
    log_value(NORMAL, "Confidence level", m_confidence);
    log_value(NORMAL, "Corresponding scaling", scale);

    // Compute minor and major axes of error ellipse
    double major = scale*std::sqrt(lambda1);
    double minor = scale*std::sqrt(lambda2);

    // Write header into logger
    log_header1(TERSE, "Generate butterfly");

    // Initialise result arrays
    m_energies.assign(m_ebounds.size(), 0.0);
    m_intensities.assign(m_ebounds.size(), 0.0);
    m_min_intensities.assign(m_ebounds.size(), 1e30);
    m_max_intensities.assign(m_ebounds.size(), 0.0);

    // Walk around the error ellipse
    int    nt = 360;
    double dt = gammalib::twopi/double(nt);
    double t  = 0.0;
    for (int i = 0; i < nt; ++i, t += dt) {

        // Pre-compute sin and cosine
        double sin_t = std::sin(t);
        double cos_t = std::cos(t);

        // Compute prefactor and index
        double prefactor = (major * cos_t * vector1[0] +
                            minor * sin_t * vector2[0]) * prefactor_scale +
                           prefactor_mean;
        double index     = (major * cos_t * vector1[1] +
                            minor * sin_t * vector2[1]) * index_scale +
                           index_mean;

        // Write prefactor and index into logger
        log_value(EXPLICIT, "Angle "+gammalib::str(t*gammalib::rad2deg),
                  gammalib::str(prefactor)+" ph/cm2/s/MeV; "+
                  gammalib::str(index));

        // Setup power law model
        GModelSpectralPlaw plaw(prefactor, index, pivot);

        // Evaluate power law at each energy and extract minimum and
        // maximum intensity
        for (int k = 0; k < m_ebounds.size(); ++k) {

            // Get the energy of current bin
            GTime   time;
            GEnergy energy = m_ebounds.elogmean(k);

            // Evaluate power law
            double intensity = plaw.eval(energy, time);

            // Get extremes
            if (intensity < m_min_intensities[k]) {
                m_min_intensities[k] = intensity;
            }
            if (intensity > m_max_intensities[k]) {
                m_max_intensities[k] = intensity;
            }

        } // endfor: looped over all energies

    } // endfor: looped over all angles

    // Setup power law model for mean prefactor and index
    GModelSpectralPlaw plaw(prefactor_mean, index_mean, pivot);

    // Compute energy and intensity vectors for mean power law
    for (int k = 0; k < m_ebounds.size(); ++k) {

        // Get the energy of current bin
        GTime   time;
        GEnergy energy = m_ebounds.elogmean(k);

        // Evaluate power law
        double intensity = plaw.eval(energy, time);

        // Store results
        m_energies[k]    = energy.MeV();
        m_intensities[k] = intensity;

    } // endfor: looped over energies

    // Write header into logger
    log_header1(TERSE, "Results");

    // Show results
    for (int k = 0; k < m_ebounds.size(); ++k) {

        // Write results into logger
        log_value(NORMAL, "Intensity at "+gammalib::str(m_energies[k])+" MeV",
                  gammalib::str(m_intensities[k])+" ("+
                  gammalib::str(m_min_intensities[k])+", "+
                  gammalib::str(m_max_intensities[k])+") ph/cm2/s/MeV");

    } // endfor: looped over energies

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
 * @brief Save butterfly diagram
 *
 * Saves the butterfly diagram into an ascii file using a column separated
 * value (CSV) format with blanks as separators.
 ***************************************************************************/
void ctbutterfly::save(void)
{
    // Write header into logger
    log_header1(TERSE, "Save Butterfly diagram");

    // Get output filename
    m_outfile = (*this)["outfile"].filename();

    // Save only if filename is non-empty
    if (!m_outfile.is_empty()) {

        // Log butterfly diagram file name
        log_value(NORMAL, "Butterfly file", m_outfile.url());

        // Create CSV table with 4 columns
        GCsv table(m_energies.size(), 4);

        // Fill CSV table
        for (int i = 0; i < m_energies.size(); ++i) {
            table.real(i, 0, m_energies[i]);
            table.real(i, 1, m_intensities[i]);
            table.real(i, 2, m_min_intensities[i]);
            table.real(i, 3, m_max_intensities[i]);
        }

        // Save CSV table
        table.save(m_outfile.url(), " ", clobber());

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
void ctbutterfly::init_members(void)
{
    // Initialise members
    m_srcname.clear();
    m_confidence  = 0.68;
    m_max_iter    = 50;
    m_apply_edisp = false;
    m_fit         = false;
    m_chatter     = static_cast<GChatter>(2);
    m_ebounds.clear();
    m_outfile.clear();

    // Initialise protected members
    m_covariance.clear();
    m_energies.clear();
    m_intensities.clear();
    m_min_intensities.clear();
    m_max_intensities.clear();

    // Set optimizer parameters
    m_opt.max_iter(m_max_iter);
    m_opt.max_stalls(10);

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
    m_srcname     = app.m_srcname;
    m_confidence  = app.m_confidence;
    m_max_iter    = app.m_max_iter;
    m_apply_edisp = app.m_apply_edisp;
    m_fit         = app.m_fit;
    m_chatter     = app.m_chatter;
    m_ebounds     = app.m_ebounds;
    m_outfile     = app.m_outfile;
 
    // Copy protected members
    m_covariance      = app.m_covariance;
    m_energies        = app.m_energies;
    m_intensities     = app.m_intensities;
    m_min_intensities = app.m_min_intensities;
    m_max_intensities = app.m_max_intensities;

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
 * @exception GException::feature_not_implemented
 *            Loading of covariance matrix is not yet implemented.
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. Most parameters are only required if no observation exists so
 * far in the observation container. In this case, a single CTA observation
 * will be added to the container, using the definition provided in the
 * parameter file.
 ***************************************************************************/
void ctbutterfly::get_parameters(void)
{
    // Setup observations from "inobs" parameter
    setup_observations(m_obs);
    
    // Setup models from "inmodel" parameter
    setup_models(m_obs, (*this)["srcname"].string());

    // Get name of test source
    m_srcname = (*this)["srcname"].string();

    // Check model name and type
    check_model();

    // Create energy boundaries from user parameters
    m_ebounds = create_ebounds();

    // Get matrix file name and load if possible
    std::string matrixfilename = (*this)["matrix"].filename();
    if (matrixfilename != "NONE") {
        std::string msg = "Loading of matrix from file not implemented yet. "
                          "Use filename = \"NONE\" to induce a recomputation of "
                          "the matrix internally.";
        throw GException::feature_not_implemented(G_GET_PARAMETERS, msg);
        // m_covariance.load(matrixfilename);
    }

    // Get remaining parameters
    m_apply_edisp = (*this)["edisp"].boolean();
    m_confidence  = (*this)["confidence"].real();
    m_fit         = (*this)["fit"].boolean();
    m_chatter     = static_cast<GChatter>((*this)["chatter"].integer());

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outfile = (*this)["outfile"].filename();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if sky model is valid
 *
 * @exception GException::invalid_value
 *            Did not find a valid model
 *
 * Checks whether the specified sky model exists and whether it is a power
 * law.
 ***************************************************************************/
void ctbutterfly::check_model(void)
{
    // Continue only if source model exists
    if (m_obs.models().contains(m_srcname)) {

        // Get relevant model and parameter for upper limit computation.
        GModels&   models = const_cast<GModels&>(m_obs.models());
        GModelSky* model  = dynamic_cast<GModelSky*>(models[m_srcname]);
        if (model == NULL) {
            std::string msg = "Source \""+m_srcname+"\" is not a sky model. "
                              "Please specify the name of a sky model for "
                              "butterfly computation.";
            throw GException::invalid_value(G_CHECK_MODEL, msg);
        }

        // Check that the spectral model is a power law
        if (model->spectral()->type() != "PowerLaw") {
            std::string msg = "\""+model->spectral()->type()+"\" cannot be "
                              "used as spectral model for an butterfly "
                              "computation. Please specify a power law model.";
            throw GException::invalid_value(G_CHECK_MODEL, msg);
        }

    } // endif: source model existed

    // Otherwise throw an exception
    else {
        std::string msg = "Source \""+m_srcname+"\" not found in model "
                          "container. Please add a source with that name "
                          "or check for possible typos.";
        throw GException::invalid_value(G_CHECK_MODEL, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute normalized eigenvectors and eigenvalues
 *
 * @param[in] a Covariance matrix element (0,0)
 * @param[in] b Covariance matrix element (0,1)
 * @param[in] c Covariance matrix element (1,0)
 * @param[in] d Covariance matrix element (1,1)
 * @param[out] lambda1 First eigenavlue
 * @param[out] lambda2 Second eigenavlue
 * @param[out] vector1 First eigenvector (normalized)
 * @param[out] vector2 Second eigenvector (normalized)
 *
 * Computes the eigenvalues and eigenvectors for a 2x2 matrix. See
 * http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
 ***************************************************************************/
void ctbutterfly::eigenvectors(const double& a,
                               const double& b,
                               const double& c,
                               const double& d,
                               double*       lambda1,
                               double*       lambda2,
                               GVector*      vector1,
                               GVector*      vector2)
{
    // Compute trace and determinant
    double trace = a + d;
    double det   = a*d - b*c;

    // Compute eigenvalues
    double arg1 = 0.5*trace;
    double arg2 = std::sqrt(0.25*trace*trace - det);
    *lambda1    = arg1 + arg2;
    *lambda2    = arg1 - arg2;

    // Compute eigenvectors
    if (c != 0.0) {
        (*vector1)[0] = *lambda1 - det;
        (*vector1)[1] = c;
        (*vector2)[0] = *lambda2 - det;
        (*vector2)[1] = c;
    }
    else if (b != 0.0) {
        (*vector1)[0] = b;
        (*vector1)[1] = *lambda1 - a;
        (*vector2)[0] = b;
        (*vector2)[1] = *lambda2 - a;
    }
    else {
        (*vector1)[0] = 1.0;
        (*vector1)[1] = 0.0;
        (*vector2)[0] = 0.0;
        (*vector2)[1] = 1.0;
    }

    // Normalize eigenvectors
    *vector1 /= norm(*vector1);
    *vector2 /= norm(*vector2);
    
    // Return
    return;
}
