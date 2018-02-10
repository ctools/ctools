/***************************************************************************
 *                       ctskymap - Sky mapping tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
#define G_GET_PARAMETERS                          "ctskymap::get_parameter()"
#define G_MAP_EVENTS                 "ctskymap::map_events(GCTAObservation*)"
#define G_MAP_BACKGROUND_IRF "ctskymap::map_background_irf(GCTAObservation*)"
#define G_MAP_BACKGROUND_RING                "ctskymap::map_background_ring("\
                                                          "GCTAObservation*)"

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
ctskymap::ctskymap(void) : ctobservation(CTSKYMAP_NAME, VERSION)
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
          ctobservation(CTSKYMAP_NAME, VERSION, obs)
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
          ctobservation(CTSKYMAP_NAME, VERSION, argc, argv)
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

    // If background subtraction is selected then compute significance map
    if (m_bkgsubtract != "NONE") {

        // Compute the significance map
        map_significance();

        // Subtract background map from counts map
        m_skymap -= m_bkgmap;
    }

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
        GFitsHDU* hdu = m_skymap.write(fits);

        // Write keywords into sky map extension
        write_ogip_keywords(hdu);
        write_hdu_keywords(hdu);

        // If background subtraction is requested then write background map
        // and significance map to FITS file
        if (m_bkgsubtract != "NONE") {
        
            // Write background map into FITS file
            hdu = m_bkgmap.write(fits);

            // Set background map extension name
            if (hdu != NULL) {
                hdu->extname("BACKGROUND");
            }

            // Write keywords into background extension
            write_ogip_keywords(hdu);
            write_hdu_keywords(hdu);

            // Write significance map into FITS file
            hdu = m_sigmap.write(fits);

            // Set significance map extension name
            if (hdu != NULL) {
                hdu->extname("SIGNIFICANCE");
            }

            // Write keywords into significance extension
            write_ogip_keywords(hdu);
            write_hdu_keywords(hdu);

        } // endif: background subtraction was requested

        // Save FITS file to disk
        fits.saveto(m_outmap, clobber());

    } // endif: filename and map were not empty

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
    m_exclmap.clear();
    m_alphamap.clear();
    m_onmap.clear();
    m_emin        = 0.0;
    m_emax        = 0.0;
    m_bkgsubtract = "NONE";
    m_roiradius   = 0.0;
    m_inradius    = 0.0;
    m_outradius   = 0.0;
    m_publish     = false;
    m_chatter     = static_cast<GChatter>(2);

    // Initialise cache
    m_solidangle.clear();
    m_dirs.clear();

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
    m_exclmap     = app.m_exclmap;
    m_alphamap    = app.m_alphamap;
    m_onmap       = app.m_onmap;
    m_emin        = app.m_emin;
    m_emax        = app.m_emax;
    m_bkgsubtract = app.m_bkgsubtract;
    m_roiradius   = app.m_roiradius;
    m_inradius    = app.m_inradius;
    m_outradius   = app.m_outradius;
    m_publish     = app.m_publish;
    m_chatter     = app.m_chatter;

    // Copy cache
    m_solidangle = app.m_solidangle;
    m_dirs       = app.m_dirs;

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

    // Get RING background parameters
    if (m_bkgsubtract == "RING") {

        // Get parameters
        m_roiradius = (*this)["roiradius"].real();
        m_inradius  = (*this)["inradius"].real();
        m_outradius = (*this)["outradius"].real();

        // Make sure that (roiradius < inradius < outradius)
        if (m_roiradius > m_inradius) {
            std::string msg("'roiradius' must be smaller than 'inradius'");
            throw GException::invalid_value(G_GET_PARAMETERS, msg);
        }
        else if (m_inradius > m_outradius) {
            std::string msg("'inradius' must be smaller than 'outradius'");
            throw GException::invalid_value(G_GET_PARAMETERS, msg);
        }

    } // endif: read RING background parameters

    // If IRF background subtraction is requested then make sure that the
    // CTA observations in the observation container have response information
    if (m_bkgsubtract != "NONE") {
        set_response(m_obs);
    }

    // Get remaining parameters
    m_publish = (*this)["publish"].boolean();
    m_chatter = static_cast<GChatter>((*this)["chatter"].integer());

    // Read ahead parameters
    if (read_ahead()) {
        m_outmap    = (*this)["outmap"].filename();
    }

    // Create background map and significance map if background subtraction
    // is requested
    if (m_bkgsubtract != "NONE") {

        // Create backgrond and significance maps
        m_bkgmap = create_map(m_obs);
        m_sigmap = create_map(m_obs);
        
        // Setup the exclusions map
        map_exclusions((*this)["inexclusion"].filename());

        // Cache the pixel solid angles and sky directions
        m_solidangle.reserve(m_bkgmap.npix());
        m_dirs.reserve(m_bkgmap.npix());
        for (int i = 0; i < m_bkgmap.npix(); ++i) {
            m_solidangle.push_back(m_bkgmap.solidangle(i));
            m_dirs.push_back(m_bkgmap.inx2dir(i));
        }

        // If doing a ring background subtraction, generate an alpha map and
        // make sure to cache the counts for each observation
        if (m_bkgsubtract == "RING") {
            m_alphamap = create_map(m_obs);
            m_onmap    = create_map(m_obs);
        }

    } // endif: background subtraction selected

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
        throw GException::no_list(G_MAP_EVENTS);
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
            num_outside_roi++;
            continue;
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
    log_header1(EXPLICIT, "Sky map");
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
 *     RING - Ring background estimation
 ***************************************************************************/
void ctskymap::map_background(GCTAObservation* obs)
{
    // Dispatch to appropriate background estimation method
    if (m_bkgsubtract == "IRF") {
        map_background_irf(obs);
    }
    else if (m_bkgsubtract == "RING") {
        map_background_ring(obs);
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
        throw GException::invalid_value(G_MAP_BACKGROUND_IRF, msg);
    }

    // Get IRF background template
    const GCTABackground* bkg = rsp->background();

    // Throw an exception if observation has no IRF background template
    if (bkg == NULL) {
        std::string msg = "No IRF background template found in instrument "
                          "response function for "+
                          get_obs_header(obs)+". Please specify an instrument "
                          "response function containing a background template.";
        throw GException::invalid_value(G_MAP_BACKGROUND_IRF, msg);
    }

    // Set minimum and maximum energy as GEnergy instances
    GEnergy emin(m_emin, "TeV");
    GEnergy emax(m_emax, "TeV");

    // Extract region of interest from observation
    GCTARoi roi = obs->roi();

    // Initialise statistics
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

        // Compute background value
        double value = bkg->rate_ebin(instdir, emin, emax);

        // Multiply background rate with livetime and solid angle
        value *= obs->livetime() * m_solidangle[i];

        // Add number of background events to background map
        m_bkgmap(i,0) += value;

        // Update total number of background events
        total += value;

    } // endfor: looped over background map pixels

    // Log background subtraction results
    log_value(NORMAL, "Events in background", int(total+0.5));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Estimates the background in sky map based on the ring background
 *        method
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No response information available for observation.
 *            No background template available in instrument response function.
 *
 * Estimates the background in the sky map by summing the events within a
 * ring centered at a given pixel's position. The pixels in the ring are
 * weighted also by the background IRF value.
 ***************************************************************************/
void ctskymap::map_background_ring(GCTAObservation* obs)
{
    // Get IRF response (to scale background counts)
    const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(obs->response());

    // Throw an exception if observation has no instrument response function
    if (rsp == NULL) {
        std::string msg = "No response information available for "+
                          get_obs_header(obs)+" to compute IRF background. "
                          "Please specify response information or use "
                          "another background subtraction method.";
        throw GException::invalid_value(G_MAP_BACKGROUND_RING, msg);
    }

    // Get IRF background template
    const GCTABackground* bkg = rsp->background();

    // Throw an exception if observation has no IRF background template
    if (bkg == NULL) {
        std::string msg = "No IRF background template found in instrument "
                          "response function for "+
                          get_obs_header(obs)+". Please specify an instrument "
                          "response function containing a background template.";
        throw GException::invalid_value(G_MAP_BACKGROUND_RING, msg);
    }

    // Set minimum and maximum energy as GEnergy instances
    GEnergy emin(m_emin, "TeV");
    GEnergy emax(m_emax, "TeV");

    // Extract region of interest from observation
    GCTARoi roi = obs->roi();

    // Extract exposure
    double exposure = obs->livetime();

    // Loop over all map pixels
    for (int i = 0; i < m_bkgmap.npix(); ++i) {

        // Get sky direction of pixel
        GSkyDir& skydir = m_dirs[i];

        // Convert sky direction in instrument direction
        GCTAInstDir instdir = obs->pointing().instdir(skydir);
    
        // If RoI is valid and instrument direction is not within RoI then
        // skip pixel
        if (roi.is_valid() && !roi.contains(instdir)) {
            continue;
        }

        // Compute background value
        m_alphamap(i,0) += bkg->rate_ebin(instdir, emin, emax) *
                           m_solidangle[i] * exposure;

        // Store the counts in this bin
        m_onmap(i,0) += m_skymap(i,0);

    } // endfor: Loop for caching bkg IRF sensitivity

    // Zero out the counts map to prepare it for the next observation
    m_skymap = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write keywords in FITS HDU
 *
 * @param[in,out] hdu Pointer to FITS HDU.
 *
 * Writes keywords in FITS HDU.
 ***************************************************************************/
void ctskymap::write_hdu_keywords(GFitsHDU* hdu) const
{
    // Continue only if pointer is valid
    if (hdu != NULL) {

        // Set keywords
        hdu->card("BKGSUB", m_bkgsubtract, "Background substraction method");
        hdu->card("E_MIN",  m_emin, "[TeV] Lower energy boundary");
        hdu->card("E_MAX",  m_emax, "[TeV] Upper energy boundary");
        hdu->card("EUNIT",  "TeV",  "Units for E_MIN and E_MAX");

    } // endif: pointer was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Computes the bin-by-bin significance if background subtraction
 *        requested
 *
 * Method computes the bin-by bin significance if a background subtraction 
 * method is specified. If method is "IRF", Poisson statistics (in the
 * Gaussian limit) are assumed. If method is "RING" a Li & Ma significance
 * is computed for each bin (as well as a map of the computed alpha values).
 ***************************************************************************/
void ctskymap::map_significance(void)
{
    // Compute significance from "RING" method (Li & Ma eq. 17)
    if (m_bkgsubtract == "RING") {

        // Log message about what is being done
        log_header1(NORMAL, "Computing Ring background map");
        log_value(NORMAL, "Total pixels to process", m_onmap.npix());

        // Store the number of bins with an inappropriate alpha parameter
        int num_bad_alpha = 0;

        // Loop through each bin in the on-counts map
        for (int i = 0; i < m_onmap.npix(); ++i) {

            // Initialise the on/off counts and alpha for this bin
            double n_on  = 0.0;
            double n_off = 0.0;
            double alpha = 0.0;

            // Since this can take a long time, keep the user updated on
            // the progress when another 10% of pixels is processed
            if (i % (m_onmap.npix()/10) == 0) {
                log_value(NORMAL, "Pixels remaining", m_onmap.npix()-i);
            }

            // Get bin coordinates
            GSkyDir& skydir = m_dirs[i];

            // Compute the alpha and counts for this bin
            compute_ring_values(m_onmap, m_alphamap, skydir,
                                n_on, n_off, alpha);

            // If alpha is zero then increment the bad alpha counter
            if (alpha == 0.0) {
                num_bad_alpha++;
            }

            // ... otherwise store the results
            else {

                // Store the on-counts and alpha-weighted off-counts
                m_skymap(i,0) = n_on;
                m_bkgmap(i,0) = alpha * n_off;

                // Compute and store significance (Li & Ma eq. 17)
                if (n_on == 0.0) {
                    m_sigmap(i,0) = -2.0 * n_off * std::log(1.0+alpha);
                }
                else {
                    m_sigmap(i,0) = (n_on < (alpha*n_off) ? -2.0 : 2.0) *
                                    (n_on  * std::log((1.0+alpha) *
                                     n_on / (alpha * (n_on+n_off))) +
                                     n_off * std::log((1.0+alpha) *
                                     n_off / (n_on+n_off)));
                }

            } // endelse: alpha was non-zero

        } // endfor: looped over all pixels

        // Log the number of bad-alpha bins
        log_value(NORMAL, "Bins with alpha=0", num_bad_alpha);

        // Now take the square root. Since some bins can be negative, first
        // take the absolute value of the map, square root that, then multiply
        // each bin by its sign to preserve the +/- significance.
        m_sigmap = sign(m_sigmap) * sqrt(abs(m_sigmap));

    } // endif: ring background selected

    // Compute significance using Poisson statistics in the Gaussian limit
    else {
        m_sigmap = (m_skymap - m_bkgmap) / sqrt(m_skymap);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Computes Non, Noff and alpha for a counts map and sensitivity map
 *
 * @param[in]  counts           Counts map
 * @param[in]  sensitivity      Sensitivity map
 * @param[in]  position         Position of the ring & roi centers
 * @param[out] non              Returned estimate of ON counts
 * @param[out] noff             Returned estimate of OFF counts
 * @param[out] alpha            Returned estimate of alpha
 * 
 * This method computes the Non counts Noff values at a given position. It 
 * also computes the alpha parameter from the passed sensitivity map.
 ***************************************************************************/
void ctskymap::compute_ring_values(const GSkyMap& counts, 
                                   const GSkyMap& sensitivity,
                                   const GSkyDir& position,
                                   double&        non,
                                   double&        noff,
                                   double&        alpha)
{
    // Reset Non and Noff
    non  = 0.0;
    noff = 0.0;

    // Initialise On and Off alpha
    double alpha_on  = 0.0;
    double alpha_off = 0.0;

    // Compute cosine of radii. We use the cosine of the radius here to
    // mimimize as much as possible the trigonometric computations. Note that
    // the cosine of an angle is maximal for an angle of zero. This explains
    // later why we chose ">=" to test whether an angle is smaller than a
    // given radius.
    double cos_roiradius = std::cos(m_roiradius * gammalib::deg2rad);
    double cos_inradius  = std::cos(m_inradius  * gammalib::deg2rad);
    double cos_outradius = std::cos(m_outradius * gammalib::deg2rad);

    // Loop over every pixel in the observation to compute Non, Noff
    for (int j = 0; j < counts.npix(); ++j) {

        // Get the index and sky direction of this pixel
        GSkyDir& skydir = m_dirs[j];

        // Only consider pixels within the outer radius of the background region
        if (position.cos_dist(skydir) >= cos_outradius) {

            // Check if pixel is inside the background region
            if ((m_exclmap(j,0) == 0.0) && position.cos_dist(skydir) < cos_inradius) {

                // Update n_off
                noff += counts(j,0);

                // Update alpha_off
                alpha_off += sensitivity(j,0);

            }

            // ... otherwise check if pixel is inside source region
            else if (position.cos_dist(skydir) >= cos_roiradius) {

                // Update n_on for significance computation
                non += counts(j,0);

                // Update alpha_on
                alpha_on += sensitivity(j,0);

            } // endif: source and background region check

        } // endif: pixel is within outer radius of background region

    } // endfor: looped over pixels

    // Compute alpha. If the off region does not have any sensitivity then
    // set Non = Noff = 0
    if (alpha_off == 0.0) {
        alpha = 0.0;
        non   = 0.0;
        noff  = 0.0;
    }
    else {
        alpha = alpha_on / alpha_off;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generates map of pixel exclusions
 *
 * @param[in] filename Exclusion file name.
 *
 * Generates a sky map of the pixels that are to be excluded from the
 * background estimation. Pixels with values different from 0 will be
 * excluded.
 ***************************************************************************/
void ctskymap::map_exclusions(const GFilename& filename)
{
    // Create exlusion map
    m_exclmap = create_map(m_obs);

    // Set all pixels to 0 (no pixel excluded)
    m_exclmap = 0.0;

    // Make sure the exclusions filename is valid
    if (is_valid_filename(filename)) {

        // Fill the exclusions based on the regions supplied
        if (filename.is_fits()) {
            map_exclusions_fits(filename);
        }
        else {
            map_exclusions_reg(filename);
        }

    } // endif: filename was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fills exclusions map from FITS image
 *
 * @param[in] filename FITS image file name.
 *
 * Sets all exclusion map pixels to 1 that correspond to non-zero pixels in
 * the exclusion sky map FITS file.
 ***************************************************************************/
void ctskymap::map_exclusions_fits(const GFilename& filename)
{
    // Load the fits image
    GSkyMap inmap(filename);

    // Loop through the individual pixels in the exclusion map
    for (int i = 0; i < m_exclmap.npix(); ++i) {

        // Get the pixel direction
        GSkyDir dir = m_exclmap.pix2dir(i);

        // Check this sky position in the fits map
        if (inmap.contains(dir) && (inmap(inmap.dir2inx(dir)) != 0.0)) {

            // Set the pixel to 1
            m_exclmap(i) = 1.0;

        } // endif: pixel,region overlap check

    } // endfor: looped over exclusion map pixels

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fills exclusions map from DS9 region file
 *
 * @param[in] filename DS9 region file name.
 *
 * Sets all exclusion map pixels to 1 that are contained in any of the DS9
 * regions.
 ***************************************************************************/
void ctskymap::map_exclusions_reg(const GFilename& filename)
{
    // Load the exclusion regions
    GSkyRegions regions(filename);

    // Loop through the individual pixels in the exclusion map
    for (int i = 0; i < m_exclmap.npix(); ++i) {

        // Get the pixel position
        GSkyDir dir = m_exclmap.pix2dir(i);

        // If pixel position overlaps with the regions
        if (regions.contains(dir)) {
            m_exclmap(i) = 1.0;
        }

    } // endfor: looped over exclusion map pixels

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
