/***************************************************************************
 *                       ctskymap - Sky mapping tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2017 by Juergen Knoedlseder                         *
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
 * @file ctskymap.cpp
 * @brief Sky mapping tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctskymap.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_INIT_MAP                 "ctskymap::init_map(GCTAObservation* obs)"
#define G_BIN_EVENTS                 "ctskymap::bin_events(GCTAObservation*)"
#define G_BKG_SUBTRACT_IRF     "ctskymap::bkg_subtract_irf(GCTAObservation*)"

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
 * Constructs an empty sky mapping tool.
 ***************************************************************************/
ctskymap::ctskymap(void) : ctobservation(CTSKYMAP_NAME, CTSKYMAP_VERSION)
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
 * Constructs sky mapping tool from an observation container.
 ***************************************************************************/
ctskymap::ctskymap(const GObservations& obs) :
          ctobservation(CTSKYMAP_NAME, CTSKYMAP_VERSION, obs)
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
 * Constructs sky mapping tool using command line arguments for user
 * parameter setting.
 ***************************************************************************/
ctskymap::ctskymap(int argc, char *argv[]) : 
          ctobservation(CTSKYMAP_NAME, CTSKYMAP_VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Sky mapping tool.
 ***************************************************************************/
ctskymap::ctskymap(const ctskymap& app) : ctobservation(app)
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
ctskymap::~ctskymap(void)
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
 * @param[in] app Sky mapping tool.
 * @return Sky mapping tool.
 ***************************************************************************/
ctskymap& ctskymap::operator=(const ctskymap& app)
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
 * @brief Clear sky mapping tool
 *
 * Clears sky mapping tool.
 ***************************************************************************/
void ctskymap::clear(void)
{
    // Free members
    free_members();
    this->ctool::free_members();
    this->ctobservation::free_members();

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
 * @brief Run the sky mapping tool
 *
 * Generates a sky map from event list by looping over all unbinned CTA
 * observation in the observation container and filling all events into
 * a sky map.
 ***************************************************************************/
void ctskymap::run(void)
{
    // Switch screen logging on in debug mode
    if (logDebug()) {
        log.cout(true);
    }

    // Get parameters
    get_parameters();

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write header into logger
    log_header1(TERSE, gammalib::number("Map observation", m_obs.size()));

    // Loop over all unbinned CTA observations in the container
    for (GCTAObservation* obs = first_unbinned_observation(); obs != NULL;
         obs = next_unbinned_observation()) {

        // Map events into sky map
        map_events(obs);

        // Compute background sky map
        map_background(obs);

        // Dispose events to free memory
        obs->dispose_events();

    } // endfor: looped over observations

    // If background subtract was required then compute significance map
    // and subtract background map from counts map
    if (m_bkgsubtract != "NONE") {

        // Compute significance map
        m_sigmap = (m_skymap - m_bkgmap) / sqrt(m_skymap);

        // Subtract background map from counts map
        m_skymap -= m_bkgmap;

    } // endif: background subtraction was requested

    // Optionally publish sky map
    if (m_publish) {
        publish();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save sky map
 *
 * Saves the sky map into a FITS file. The FITS file name is specified by the
 * @p outname parameter.
 ***************************************************************************/
void ctskymap::save(void)
{
    // Write header
    log_header1(TERSE, "Save sky map");

    // Get sky map filename
    m_outmap  = (*this)["outmap"].filename();

    // Save sky map if filename and the map are not empty
    if (!m_outmap.is_empty() && !m_skymap.is_empty()) {

        // Create empty FITS file
        GFits fits;

        // Write sky map into FITS file
        m_skymap.write(fits);

        // If background subtraction is requested then write background map
        // and significance map to FITS file
        if (m_bkgsubtract != "NONE") {
        
            // Write background map into FITS file
            m_bkgmap.write(fits);

            // Set background map extension name
            fits[fits.size()-1]->extname("BACKGROUND");

            // Write significance map into FITS file
            m_sigmap.write(fits);

            // Set significance map extension name
            fits[fits.size()-1]->extname("SIGNIFICANCE");

        } // endif: background subtraction was requested

        // Save FITS file to disk
        fits.saveto(m_outmap, clobber());

    }

    // Write into logger what has been done
    std::string fname = (m_outmap.is_empty()) ? "NONE" : m_outmap.url();
    if (m_skymap.is_empty()) {
        fname.append(" (map is empty, no file created)");
    }
    log_value(NORMAL, "Sky map file", fname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Publish sky map
 *
 * @param[in] name Sky map name.
 ***************************************************************************/
void ctskymap::publish(const std::string& name)
{
    // Write header into logger
    log_header1(TERSE, "Publish sky map");

    // Set default name if user name is empty
    std::string user_name(name);
    if (user_name.empty()) {
        user_name = CTSKYMAP_NAME;
    }

    // Write sky map name into logger
    log_value(NORMAL, "Sky map name", user_name);

    // Publish sky map
    m_skymap.publish(user_name);

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
void ctskymap::init_members(void)
{
    // Initialise members
    m_skymap.clear();
    m_bkgmap.clear();
    m_sigmap.clear();
    m_emin        = 0.0;
    m_emax        = 0.0;
    m_bkgsubtract = "NONE";
    m_publish     = false;
    m_chatter     = static_cast<GChatter>(2);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctskymap::copy_members(const ctskymap& app)
{
    // Copy attributes
    m_skymap      = app.m_skymap;
    m_bkgmap      = app.m_bkgmap;
    m_sigmap      = app.m_sigmap;
    m_emin        = app.m_emin;
    m_emax        = app.m_emax;
    m_bkgsubtract = app.m_bkgsubtract;
    m_publish     = app.m_publish;
    m_chatter     = app.m_chatter;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctskymap::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. Most parameters are only required if no observation exists so
 * far in the observation container. In this case, a single CTA observation
 * will be added to the container, using the definition provided in the
 * parameter file.
 ***************************************************************************/
void ctskymap::get_parameters(void)
{
    // Setup observations from "inobs" parameter. Do not request response
    // information and do not accept counts cubes.
    setup_observations(m_obs, false, true, false);

    // Create sky map based on task parameters
    m_skymap = create_map(m_obs);

    // Get further parameters
    m_emin        = (*this)["emin"].real();
    m_emax        = (*this)["emax"].real();
    m_bkgsubtract = (*this)["bkgsubtract"].string();

    // If IRF background subtraction is requested then make sure that the
    // CTA observations in the observation container have response information
    if (m_bkgsubtract == "IRF") {
        set_response(m_obs);
    }

    // Get remaining parameters
    m_publish     = (*this)["publish"].boolean();
    m_chatter     = static_cast<GChatter>((*this)["chatter"].integer());

    // Read ahead parameters
    if (read_ahead()) {
        m_outmap  = (*this)["outmap"].filename();
    }

    // Create background map and significance map if background subtraction
    // is requested
    if (m_bkgsubtract != "NONE") {
        m_bkgmap = create_map(m_obs);
        m_sigmap = create_map(m_obs);
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Map events into a sky map
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::no_list
 *            No event list found in observation.
 *
 * This method maps the events found in a CTA events list into a sky map.
 ***************************************************************************/
void ctskymap::map_events(GCTAObservation* obs)
{
    // Get non-const pointer on a CTA event list
    GCTAEventList* events = dynamic_cast<GCTAEventList*>
                            (const_cast<GEvents*>(obs->events()));

    // Make sure that the observation holds a CTA event list. If this
    // is not the case then throw an exception.
    if (events == NULL) {
        throw GException::no_list(G_BIN_EVENTS);
    }

    // Setup energy range covered by data
    GEnergy  emin;
    GEnergy  emax;
    GEbounds ebds;
    emin.TeV(m_emin);
    emax.TeV(m_emax);

    // Initialise binning statistics
    int num_outside_roi    = 0;
    int num_outside_map    = 0;
    int num_outside_erange = 0;
    int num_in_map         = 0;

    // Extract region of interest from observation
    GCTARoi roi = obs->roi();

    // Fill sky map
    for (int i = 0; i < events->size(); ++i) {

        // Get event
        GCTAEventAtom* event = (*events)[i];

        // Skip if energy is out of range
        if (event->energy() < emin || event->energy() > emax) {
            num_outside_erange++;
            continue;
        }

        // Determine sky pixel
        GCTAInstDir* inst  = (GCTAInstDir*)&(event->dir());
        GSkyDir      dir   = inst->dir();
        GSkyPixel    pixel = m_skymap.dir2pix(dir);

        // Skip if pixel is out of range
        if (pixel.x() < -0.5 || pixel.x() > (m_skymap.nx() - 0.5) ||
            pixel.y() < -0.5 || pixel.y() > (m_skymap.ny() - 0.5)) {
            num_outside_map++;
            continue;
        }

        // If RoI is valid then skip if  instrument direction is not within RoI
        if (roi.is_valid() && !roi.contains(*inst)) {
            continue;
            num_outside_roi++;
        }

        // Fill event in skymap
        m_skymap(pixel, 0) += 1.0;
        num_in_map++;

    } // endfor: looped over all events

    // Log binning results
    log_value(NORMAL, "Events in list", obs->events()->size());
    log_value(NORMAL, "Events in map", num_in_map);
    log_value(NORMAL, "Events outside RoI", num_outside_roi);
    log_value(NORMAL, "Events outside map area", num_outside_map);
    log_value(NORMAL, "Events outside energies", num_outside_erange);

    // Write sky map into header
    log_header2(EXPLICIT, "Sky map");
    log_string(EXPLICIT, m_skymap.print(m_chatter));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Estimates the background in sky map for observation
 *
 * @param[in] obs CTA observation.
 *
 * Estimates the background in the sky map for a given observation and adds
 * this estimate to the background sky map. The background estimation method
 * is specified by the "bkgsubtract" parameter which can take the following
 * values:
 *
 *     NONE - No background estimation
 *     IRF  - Background estimation based on IRF template
 ***************************************************************************/
void ctskymap::map_background(GCTAObservation* obs)
{
    // Dispatch to appropriate background estimation method
    if (m_bkgsubtract == "IRF") {
        map_background_irf(obs);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Estimates the background in sky map based on IRF template
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No response information available for observation.
 *            No background template available in instrument response function.
 *
 * Estimates the background in the sky map using the IRF template for a given
 * observation and adds this estimate to the background sky map. The IRF
 * template is integrated numerically in energy.
 ***************************************************************************/
void ctskymap::map_background_irf(GCTAObservation* obs)
{
    // Get IRF response
    const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(obs->response());

    // Throw an exception if observation has no instrument response function
    if (rsp == NULL) {
        std::string msg = "No response information available for "+
                          get_obs_header(obs)+" to compute IRF background. "
                          "Please specify response information or use "
                          "another background subtraction method.";
        throw GException::invalid_value(G_BKG_SUBTRACT_IRF, msg);
    }

    // Get IRF background template
    const GCTABackground* bkg = rsp->background();

    // Throw an exception if observation has no IRF background template
    if (bkg == NULL) {
        std::string msg = "No IRF background template found in instrument "
                          "response function for "+
                          get_obs_header(obs)+". Please specify an instrument "
                          "response function containing a background template.";
        throw GException::invalid_value(G_BKG_SUBTRACT_IRF, msg);
    }

    // Compute natural logarithm of energy range in MeV
    double lnEmin = std::log(m_emin * 1.0e6);
    double lnEmax = std::log(m_emax * 1.0e6);

    // Extract region of interest from observation
    GCTARoi roi = obs->roi();

    // Initialise statistics
    int    calls = 0;
    double total = 0.0;

    // Loop over all background map pixels
    for (int i = 0; i < m_bkgmap.npix(); ++i) {

        // Get sky direction of pixel
        GSkyDir skydir = m_bkgmap.inx2dir(i);

        // Convert sky direction in instrument direction
        GCTAInstDir instdir = obs->pointing().instdir(skydir);

        // If RoI is valid and instrument direction is not within RoI then
        // skip pixel
        if (roi.is_valid() && !roi.contains(instdir)) {
            continue;
        }

        // Setup integration function
        ctskymap::irf_kern integrand(bkg, &instdir);
        GIntegral          integral(&integrand);

        // Set precision (has been carefully adjusted using a test simulation
        // over the energy range 20 GeV - 120 TeV)
        integral.eps(1.0e-6);

        // Do Romberg integration
        double value = integral.romberg(lnEmin, lnEmax);

        // Update number of background function calls
        calls += integral.calls();
        
        // Multiply background rate with livetime and solid angle
        value *= obs->livetime() * m_bkgmap.solidangle(i);

        // Add number of background events to background map
        m_bkgmap(i,0) += value;

        // Update total number of background events
        total += value;

    } // endfor: looped over background map pixels
  
    // Log background subtraction results
    log_value(NORMAL, "Events in background", int(total+0.5));
    log_value(NORMAL, "Background evaluations", calls);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Estimates the background in sky map based on IRF template
 *
 * @param[in] lnE Natural logarithm of energy in MeV.
 ***************************************************************************/
double ctskymap::irf_kern::eval(const double& lnE)
{
    // Get log10 of energy in TeV
    double logE = lnE * gammalib::inv_ln10 - 6;

    // Get function value
    double value = (*m_bgd)(logE, m_dir->detx(), m_dir->dety());

    // Correct for variable substitution
    value *= std::exp(lnE);

    // Return value
    return value;
}

