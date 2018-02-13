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
#define G_RING_KERNEL               "ctskymap::ring_kernel(double&, double&)"

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

    // Setup maps
    setup_maps();

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write header into logger
    log_header1(TERSE, gammalib::number("Find unbinned observation",
                m_obs.size()));

    // Find all unbinned CTA observations in m_obs
    std::vector<GCTAObservation*> obs_list(0);
    for (GCTAObservation* obs = first_unbinned_observation(); obs != NULL;
            obs = next_unbinned_observation()) {

        // Push observation into list
        obs_list.push_back(obs);

        // Write message
        std::string msg = " Including unbinned "+obs->instrument()+
                          " observation";
        log_string(NORMAL, msg);
    }

    // Write header into logger
    log_header1(TERSE, gammalib::number("Map observation", m_obs.size()));

    // Loop over all unbinned CTA observations in the container
    #pragma omp parallel for
    for (int i = 0; i < obs_list.size(); ++i) {

        // Get pointer to observation
        GCTAObservation* obs = obs_list[i];

        // Map events into sky map
        map_events(obs);

        // Compute background sky map (only if m_bkgsubtract != "NONE")
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

    } // endif: background subtraction is selected

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
    // Initialise User parameters
    m_outmap.clear();
    m_inexclusion.clear();
    m_binsz       = 0.0;
    m_emin        = 0.0;
    m_emax        = 0.0;
    m_bkgsubtract = "NONE";
    m_roiradius   = 0.0;
    m_inradius    = 0.0;
    m_outradius   = 0.0;
    m_usefft      = true;
    m_publish     = false;
    m_chatter     = static_cast<GChatter>(2);

    // Initialise members
    m_skymap.clear();
    m_bkgmap.clear();
    m_sigmap.clear();
    m_exclmap.clear();

    // Initialise cache
    m_solidangle.clear();
    m_dirs.clear();
    m_cos_roiradius = 1.0;
    m_cos_inradius  = 1.0;
    m_cos_outradius = 1.0;

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
    // Copy User parameters
    m_outmap      = app.m_outmap;
    m_inexclusion = app.m_inexclusion;
    m_binsz       = app.m_binsz;
    m_emin        = app.m_emin;
    m_emax        = app.m_emax;
    m_bkgsubtract = app.m_bkgsubtract;
    m_roiradius   = app.m_roiradius;
    m_inradius    = app.m_inradius;
    m_outradius   = app.m_outradius;
    m_usefft      = app.m_usefft;
    m_publish     = app.m_publish;
    m_chatter     = app.m_chatter;

    // Copy members
    m_skymap  = app.m_skymap;
    m_bkgmap  = app.m_bkgmap;
    m_sigmap  = app.m_sigmap;
    m_exclmap = app.m_exclmap;

    // Copy cache
    m_solidangle    = app.m_solidangle;
    m_dirs          = app.m_dirs;
    m_cos_roiradius = app.m_cos_roiradius;
    m_cos_inradius  = app.m_cos_inradius;
    m_cos_outradius = app.m_cos_outradius;

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

    // Create sky map. Note that the method reads the relevant User parameters.
    m_skymap = create_map(m_obs);

    // Get further parameters
    m_binsz       = (*this)["binsz"].real();
    m_emin        = (*this)["emin"].real();
    m_emax        = (*this)["emax"].real();
    m_bkgsubtract = (*this)["bkgsubtract"].string();

    // Get RING background parameters
    if (m_bkgsubtract == "RING") {

        // Get parameters
        m_roiradius   = (*this)["roiradius"].real();
        m_inradius    = (*this)["inradius"].real();
        m_outradius   = (*this)["outradius"].real();
        m_inexclusion = (*this)["inexclusion"].filename();
        m_usefft      = (*this)["usefft"].boolean();

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

    // If background subtraction is requested then make sure that the CTA
    // observations in the observation container have response information
    if (m_bkgsubtract != "NONE") {
        set_response(m_obs);
    }

    // Get remaining parameters
    m_publish = (*this)["publish"].boolean();
    m_chatter = static_cast<GChatter>((*this)["chatter"].integer());

    // Read ahead parameters
    if (read_ahead()) {
        m_outmap = (*this)["outmap"].filename();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup maps
 *
 * Setup maps for sky map generation.
 ***************************************************************************/
void ctskymap::setup_maps(void)
{
    // Clear cache
    m_solidangle.clear();
    m_dirs.clear();

    // Create background map and significance map if background subtraction
    // is requested
    if (m_bkgsubtract != "NONE") {

        // Create background and significance maps
        m_bkgmap = create_map(m_obs);
        m_sigmap = create_map(m_obs);

        // Setup the exclusions map
        map_exclusions(m_inexclusion);

        // Cache the pixel solid angles and sky directions
        m_solidangle.reserve(m_bkgmap.npix());
        m_dirs.reserve(m_bkgmap.npix());
        for (int i = 0; i < m_bkgmap.npix(); ++i) {
            m_solidangle.push_back(m_bkgmap.solidangle(i));
            m_dirs.push_back(m_bkgmap.inx2dir(i));
        }

    } // endif: background subtraction selected

    // Compute cosine of ring radii for RING background method
    if (m_bkgsubtract == "RING") {
        m_cos_roiradius = std::cos(m_roiradius * gammalib::deg2rad);
        m_cos_inradius  = std::cos(m_inradius  * gammalib::deg2rad);
        m_cos_outradius = std::cos(m_outradius * gammalib::deg2rad);
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
            map_exclusions_region(filename);
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
void ctskymap::map_exclusions_region(const GFilename& filename)
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
        #pragma omp critical(ctskymap_map_events)
        m_skymap(pixel) += 1.0;
        num_in_map++;

    } // endfor: looped over all events

    // Log binning results
    #pragma omp critical(ctskymap_map_events)
    {
        log_header3(TERSE, get_obs_header(obs));
        log_value(NORMAL, "Events in list", obs->events()->size());
        log_value(NORMAL, "Events in map", num_in_map);
        log_value(NORMAL, "Events outside RoI", num_outside_roi);
        log_value(NORMAL, "Events outside map area", num_outside_map);
        log_value(NORMAL, "Events outside energies", num_outside_erange);

        // Write sky map into header
        log_header1(EXPLICIT, "Sky map");
        log_string(EXPLICIT, m_skymap.print(m_chatter));
    }

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
    const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>
                                 (obs->response());

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

    // Extract region of interest center and radius from observation and
    // compute the cosine of the RoI radius. This allows to mimimize as
    // much as possible the trigonometric computations. Note that the cosine
    // of an angle is maximal for an angle of zero. This explains later why
    // we chose ">=" to test whether an angle is smaller than a given radius.
    GSkyDir roi_centre     = obs->roi().centre().dir();
    double  cos_roi_radius = std::cos(obs->roi().radius() * gammalib::deg2rad);

    // Extract exposure
    double exposure = obs->livetime();

    // Initialise statistics
    double total = 0.0;

    // Loop over all background map pixels
    for (int i = 0; i < m_bkgmap.npix(); ++i) {

        // Get sky direction of pixel
        GSkyDir& skydir = m_dirs[i];

        // Skip pixel if it is not within the RoI
        if (roi_centre.cos_dist(skydir) < cos_roi_radius) {
            continue;
        }

        // Convert sky direction to instrument direction
        GCTAInstDir instdir = obs->pointing().instdir(skydir);

        // Compute background value
        double value = bkg->rate_ebin(instdir, emin, emax);

        // Multiply background rate with livetime and solid angle
        value *= exposure * m_solidangle[i];

        // Add number of background events to background map
        #pragma omp critical(ctskymap_map_background_irf)
        m_bkgmap(i) += value;

        // Update total number of background events
        total += value;

    } // endfor: looped over background map pixels

    // Log background subtraction results
    #pragma omp critical(ctskymap_map_events)
    {
        log_value(NORMAL, "Events in background", int(total+0.5));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute the pixel significance if background subtraction is
 *        requested
 *
 * Computes the pixel significance if a background subtraction method is
 * specified.
 *
 * For the @c IRF method Poisson statistic (in the Gaussian limit) is
 * assumed, and the significance \f$\sigma_i\f$ is computed using
 *
 * \f[\sigma_i = \frac{N_i - B_i}{\sqrt{N_i}}\f]
 *
 * where
 * \f$N_i\f$ is the number of observed counts and
 * \f$B_i\f$ is the estimated number of background counts.
 *
 * For the @c RING method the Li & Ma significance is computed for each
 * pixel, using the @c roiradius, @c inradius and @c outradius parameters
 * to specify the radius of the On region and the Off background ring,
 * respectively. Depending on the @c usefft parameter, either an FFT is used
 * to compute the number of events and the acceptance of the On and Off
 * regions, or a direct computation is performed. While the former is faster
 * it is less accurate since it assumes Euclidean distances in the sky map.
 * The latter is slower but compute exact distance between pixels of the sky
 * map. Use FFT if the sky map is close to a cartesian grid, and use the
 * direct method if the sky map shows important distortions.
 ***************************************************************************/
void ctskymap::map_significance(void)
{
    // Compute significance from "RING" method (Li & Ma eq. 17)
    if (m_bkgsubtract == "RING") {
        if (m_usefft) {
            map_significance_ring_fft();
        }
        else {
            map_significance_ring_direct();
        }
    }

    // Compute significance using Poisson statistics in the Gaussian limit
    else {
        m_sigmap = (m_skymap - m_bkgmap) / sqrt(m_skymap);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute the pixel significance for RING background
 *
 * Computes the Li & Ma significance for each sky map pixel and replaces
 * the sky and background maps by the On- and Off-count maps. The computation
 * is done by computing the exact distances between pixels in the sky map.
 ***************************************************************************/
void ctskymap::map_significance_ring_direct(void)
{
    // Log message about what is being done
    log_header1(NORMAL, "Computing ring background map (direct method)");
    log_value(NORMAL, "Total pixels to process", m_skymap.npix());

    // Initialise On- and Off-count maps
    GSkyMap onmap  = m_skymap;
    GSkyMap offmap = m_skymap;
    onmap  = 0.0;
    offmap = 0.0;

    // Initialise statistics
    int num_bad_alpha = 0;

    // Set progress logging frequency. By default the user is updated on the
    // progress when another 10% of pixels is processed, however, if there
    // are many pixels the user is updated each 100000 pixels.
    int n_logpix = m_skymap.npix() / 10;
    if (n_logpix > 100000) {
        n_logpix = 100000;
    }

    // Loop over sky map pixels
    #pragma omp parallel for
    for (int i = 0; i < m_skymap.npix(); ++i) {

        // Initialise the on/off counts and alpha for this bin
        double n_on  = 0.0;
        double n_off = 0.0;
        double alpha = 0.0;

        // Update user about progress (only if there are more than 10000 pixels)
        if (i % n_logpix == 0 && n_logpix > 10000) {
            #pragma omp critical(ctskymap_map_significance_ring)
            log_value(NORMAL, "Pixels remaining", m_skymap.npix()-i);
        }

        // Compute the alpha and counts for this bin
        compute_ring_values(i, m_skymap, m_bkgmap, n_on, n_off, alpha);

        // Store the On and alpha-weighted Off-counts
        onmap(i)  = n_on;
        offmap(i) = alpha * n_off;

        // If alpha is zero then increment the bad alpha counter
        if (alpha == 0.0) {
            num_bad_alpha++;
            #pragma omp critical(ctskymap_map_significance_ring)
            m_sigmap(i) = 0.0;
        }

        // ... otherwise store the results
        else {

            // Compute and store significance (Li & Ma eq. 17)
            double sigma = sigma_li_ma(n_on, n_off, alpha);
            #pragma omp critical(ctskymap_map_significance_ring)
            m_sigmap(i) = sigma;

        } // endelse: alpha was non-zero

    } // endfor: looped over all pixels

    // Log the number of bad-alpha bins
    log_value(NORMAL, "Pixels with alpha=0", num_bad_alpha);

    // Now take the square root. Since some bins can be negative, first
    // take the absolute value of the map, square root that, then multiply
    // each bin by its sign to preserve the +/- significance.
    m_sigmap = sign(m_sigmap) * sqrt(abs(m_sigmap));

    // Terminate by storing the on-counts map as sky map and the off-counts
    // map as background map
    m_skymap = onmap;
    m_bkgmap = offmap;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute the pixel significance for RING background
 *
 * Computes the Li & Ma significance for each sky map pixel and replaces
 * the sky and background maps by the On- and Off-count maps. The computation
 * is done using a FFT.
 ***************************************************************************/
void ctskymap::map_significance_ring_fft(void)
{
    // Log message about what is being done
    log_header1(NORMAL, "Computing ring background map (FFT method)");

    // Initialise statistics
    int num_bad_alpha = 0;

    // Set-up sky and acceptance maps with excluded pixels according to the
    // exclusion map
    GSkyMap skymap     = m_skymap;
    GSkyMap acceptance = m_bkgmap;
    for (int i = 0; i < m_skymap.npix(); ++i) {
        if (m_exclmap(i) != 0.0) {
            skymap(i)     = 0.0;
            acceptance(i) = 0.0;
        }
    }

    // Convolve maps for On-region and Off-region
    GSkyMap on_counts  = ring_convolve(skymap,     0.0,        m_roiradius);
    GSkyMap on_alpha   = ring_convolve(acceptance, 0.0,        m_roiradius);
    GSkyMap off_counts = ring_convolve(skymap,     m_inradius, m_outradius);
    GSkyMap off_alpha  = ring_convolve(acceptance, m_inradius, m_outradius);

    // Loop over sky map pixels
    for (int i = 0; i < m_skymap.npix(); ++i) {

        // Get On-counts, Off-counts and alpha value for this bin
        double n_on  = on_counts(i);
        double n_off = off_counts(i);
        double alpha;
        if (off_alpha(i) == 0.0) {
            alpha = 0.0;
            n_on  = 0.0;
            n_off = 0.0;
        }
        else {
            alpha = on_alpha(i) / off_alpha(i);
        }

        // Store the On and alpha-weighted Off-counts
        on_counts(i)  = n_on;
        off_counts(i) = alpha * n_off;

        // If alpha is zero then increment the bad alpha counter
        if (alpha == 0.0) {
            m_sigmap(i) = 0.0;
            num_bad_alpha++;
        }

        // ... otherwise compute and store the results
        else {
            m_sigmap(i) = sigma_li_ma(n_on, n_off, alpha);
        }

    } // endfor: looped over all pixels

    // Log the number of bad-alpha bins
    log_value(NORMAL, "Total number of pixels", m_skymap.npix());
    log_value(NORMAL, "Pixels with alpha=0", num_bad_alpha);

    // Now take the square root. Since some bins can be negative, first
    // take the absolute value of the map, square root that, then multiply
    // each bin by its sign to preserve the +/- significance.
    m_sigmap = sign(m_sigmap) * sqrt(abs(m_sigmap));

    // Terminate by storing the on-counts map as sky map and the off-counts
    // map as background map
    m_skymap = on_counts;
    m_bkgmap = off_counts;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Computes Non, Noff and alpha for a counts map and sensitivity map
 *
 * @param[in]  ipixel      Sky map pixel to consider.
 * @param[in]  counts      Counts map.
 * @param[in]  background  Background map.
 * @param[out] non         Returned estimate of On-counts.
 * @param[out] noff        Returned estimate of Off-counts.
 * @param[out] alpha       Returned estimate of alpha.
 * 
 * Computes the On- and Off-counts values at a given position from in counts
 * map and the alpha parameter from the background map.
 *
 * To speed-up the computations the method considers only the pixels within
 * a bounding box that comprises the outer background ring radius. The method
 * considers wrapping around of pixel indices in the Right Ascension /
 * Galactic longitude direction.
 *
 * If the alpha of the ring region is zero the values of @p non, @p noff,
 * and @p alpha will be set to zero.
 ***************************************************************************/
void ctskymap::compute_ring_values(const int&     ipixel,
                                   const GSkyMap& counts,
                                   const GSkyMap& background,
                                   double&        non,
                                   double&        noff,
                                   double&        alpha)
{
    // Initialise return values (non, noff and alpha)
    non   = 0.0;
    noff  = 0.0;
    alpha = 0.0;

    // Initialise On and Off alpha
    double alpha_on  = 0.0;
    double alpha_off = 0.0;

    // Get sky direction of pixel
    GSkyDir& position = m_dirs[ipixel];

    // Get pixel bounding box
    int ix_start;
    int ix_stop;
    int iy_start;
    int iy_stop;
    ring_bounding_box(ipixel, ix_start, ix_stop, iy_start, iy_stop);

    // Get number of x pixels in sky map
    int nx = counts.nx();

    // Loop over x pixels of bounding box
    for (int ix_nowrap = ix_start; ix_nowrap < ix_stop; ++ix_nowrap) {

        // Compute x pixel index by wrapping the pixel index into the
        // interval [0,nx[
        int ix = ix_nowrap;
        if (ix < 0) {
            ix += nx;
        }
        else if (ix >= nx) {
            ix -= nx;
        }

        // Initialise pixel index
        int i = ix + iy_start * nx;

        // Loop over y pixels of bounding box
        for (int iy = iy_start; iy < iy_stop; ++iy, i += nx) {

            // Get the index and sky direction of this pixel
            GSkyDir& skydir = m_dirs[i];

            // Only consider pixels within the outer radius of the background
            // region
            if (position.cos_dist(skydir) >= m_cos_outradius) {

                // Check if pixel is inside the background region
                if ((m_exclmap(i) == 0.0) &&
                    (position.cos_dist(skydir) < m_cos_inradius)) {

                    // Update n_off
                    noff += counts(i);

                    // Update alpha_off
                    alpha_off += background(i);

                }

                // ... otherwise check if pixel is inside source region
                else if (position.cos_dist(skydir) >= m_cos_roiradius) {

                    // Update n_on
                    non += counts(i);

                    // Update alpha_on
                    alpha_on += background(i);

                } // endif: source and background region check

            } // endif: pixel is within outer radius of background region

        } // endfor: looped over y pixels

    } // endfor: looped over x pixels

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
 * @brief Computes bounding box for RING background computation
 *
 * @param[in]  ipixel Sky map pixel to consider.
 * @param[out] ix1    Index of first pixel in x.
 * @param[out] ix2    Index after last pixel in x.
 * @param[out] iy1    Index of first pixel in y.
 * @param[out] iy2    Index after last pixel in y.
 *
 * Computes the bounding box that contained the background ring for a
 * specific pixel for the RING background method.
 *
 * The method determines the local pixel scale for the requested pixel and
 * draws a bounding box with 1.5 times the outer background ring radius
 * around the pixel. In the x direction the pixel indices are unconstrained,
 * and the client has to assure that the pixel value is comprised within
 * the validity range (this is required to handle the longitude wrap around).
 * In the y direction the pixel indices are constrained to [0,ny], where
 * ny is the number of y pixels in the sky map.
 ***************************************************************************/
void ctskymap::ring_bounding_box(const int& ipixel, int& ix1, int& ix2,
                                                    int& iy1, int& iy2) const
{
    // Get number of pixels in x and y direction
    int nx = m_skymap.nx();
    int ny = m_skymap.ny();

    // Get x and y index of pixel to consider
    int ix0 = ipixel % nx;
    int iy0 = ipixel / nx;

    // Compute pixel increment in x-direction
    double dx = 0.0;
    int    ix = (ix0 > 0) ? ix0 - 1 : ix0 + 1;
    if (ix < nx) {
        dx = m_dirs[ipixel].dist_deg(m_dirs[ix+iy0*nx]);
    }
    if (dx < m_binsz) {
        dx = m_binsz;
    }

    // Compute pixel increment in y-direction
    double dy = 0.0;
    int    iy = (iy0 > 0) ? iy0 - 1 : iy0 + 1;
    if (iy < ny) {
        dx = m_dirs[ipixel].dist_deg(m_dirs[ix0+iy*nx]);
    }
    if (dy < m_binsz) {
        dy = m_binsz;
    }

    // Compute bounding box half size in x and y. The outer radius is
    // multiplied by 1.5 to have some margin in case of map distortions
    int nx_bbox  = int(1.5 * m_outradius / dx);
    int ny_bbox  = int(1.5 * m_outradius / dy);

    // Compute index range for x. We allow that indices are smaller than 0
    // or equal or larger than nx to handle wrap around, which needs to be
    // done in the client method.
    if (2*nx_bbox < nx) {
        ix1 = ix0 - nx_bbox;
        ix2 = ix0 + nx_bbox;
    }
    else {
        ix1 = 0;
        ix2 = nx;
    }

    // Compute index range for y. The index range is restricted to [0,ny].
    if (2*ny_bbox < ny) {
        iy1 = iy0 - ny_bbox;
        iy2 = iy0 + ny_bbox;
        if (iy1 < 0) {
            iy1 = 0;
        }
        if (iy2 > ny) {
            iy2 = ny;
        }
    }
    else {
        iy1 = 0;
        iy2 = ny;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return FFT kernel for background ring
 *
 * @param[in] map  Sky map to be convolved with a background ring.
 * @param[in] rmin Minimum ring radius (degrees).
 * @param[in] rmax Maximum ring radius (degrees).
 * @return FFT kernel.
 *
 * Computes the FFT kernel for a backgrund ring.
 ***************************************************************************/
GSkyMap ctskymap::ring_convolve(const GSkyMap& map, const double& rmin,
                                                    const double& rmax) const
{
    // Copy input map
    GSkyMap conv_map = map;

    // Get FFT of ring kernel
    GFft fft_kernel(ring_kernel(rmin, rmax));

    // Extract sky map into Ndarray
    GNdarray array(conv_map.nx(), conv_map.ny());
    double *dst = array.data();
    for (int i = 0; i < conv_map.npix(); ++i) {
        *dst++ = conv_map(i);
    }

    // FFT of sky map
    GFft fft_array = GFft(array);

    // Multiply FFT of sky map with FFT of kernel
    GFft fft_smooth = fft_array * fft_kernel;

    // Backward transform sky map
    GNdarray smooth = fft_smooth.backward();

    // Insert sky map
    double *src = smooth.data();
    for (int i = 0; i < conv_map.npix(); ++i) {
        conv_map(i) = *src++;
    }

    // Return convolved map
    return conv_map;
}


/***********************************************************************//**
 * @brief Return FFT kernel for background ring
 *
 * @param[in] rmin Minimum ring radius (degrees).
 * @param[in] rmax Maximum ring radius (degrees).
 * @return FFT kernel.
 *
 * @exception GException::invalid_value
 *            Sky map is not a WCS projection.
 *
 * Computes the FFT kernel for a background ring. The computation is done
 * assuming that the sky map represents a cartesian grid, and distances are
 * computed in that grid using Euclidean distances. The @p rmin and @p rmax
 * arguments refer to a ring radius in these Euclidean distances.
 ***************************************************************************/
GNdarray ctskymap::ring_kernel(const double& rmin, const double& rmax) const
{
    // Get pointer on WCS projection
    const GWcs* wcs = dynamic_cast<const GWcs*>(m_skymap.projection());
    if (wcs == NULL) {
        std::string msg = "Sky map is not a WCS projection. Method is only "
                          "valid for WCS projections.";
        throw GException::invalid_value(G_RING_KERNEL, msg);
    }

    // Get X and Y step size
    double dx = wcs->cdelt(0);
    double dy = wcs->cdelt(1);

    // Initialise kernel
    GNdarray kern(m_skymap.nx(), m_skymap.ny());

    // Fill kernel
    for (int ix1 = 0, ix2 = m_skymap.nx(); ix1 < m_skymap.nx(); ++ix1, --ix2) {
        double x   = ix1 * dx;
        double xqs = x * x;
        for (int iy1 = 0, iy2 = m_skymap.ny(); iy1 < m_skymap.ny(); ++iy1, --iy2) {
            double y = iy1 * dy;
            double r = std::sqrt(xqs + y*y);
            if ((r >= rmin) && (r <= rmax)) {
                kern(ix1,iy1) += 1.0;
                if (ix2 < m_skymap.nx()) {
                    kern(ix2,iy1) += 1.0;
                }
                if (iy2 < m_skymap.ny()) {
                    kern(ix1,iy2) += 1.0;
                }
                if ((ix2 < m_skymap.nx()) && (iy2 < m_skymap.ny())) {
                    kern(ix2,iy2) += 1.0;
                }
            }
        }
    }

    // Return kernel
    return kern;
}


/***********************************************************************//**
 * @brief Compute significance following Li & Ma
 *
 * @param[in] n_on Number of On-counts.
 * @param[in] n_off Number of Off-counts.
 * @param[in] alpha Alpha parameters.
 * @return Significance.
 *
 * Computes the significance following Li & Ma, Equation (17).
 ***************************************************************************/
double ctskymap::sigma_li_ma(const double& n_on,
                             const double& n_off,
                             const double& alpha) const
{
    // Allocate result
    double sigma;

    // Handle special case of no On-counts
    if (n_on == 0.0) {
        sigma = -2.0 * n_off * std::log(1.0+alpha);
    }

    // ... otherwise handle general case
    else {
        sigma = (n_on < (alpha*n_off) ? -2.0 : 2.0) *
                (n_on  * std::log((1.0+alpha) * n_on / (alpha * (n_on+n_off))) +
                 n_off * std::log((1.0+alpha) * n_off / (n_on+n_off)));
    }

    // Return sigma
    return sigma;
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

        // Write RING background keywords
        if (m_bkgsubtract == "RING") {
            hdu->card("ROIRAD", m_roiradius, "[deg] Source region radius");
            hdu->card("INRAD",  m_inradius,  "[deg] Inner background ring radius");
            hdu->card("OUTRAD", m_outradius, "[deg] Outer background ring radius");
        }

    } // endif: pointer was valid

    // Return
    return;
}
