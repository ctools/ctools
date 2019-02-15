/***************************************************************************
 *                ctfindvar - Time variability search tool                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018-2019 by Simon Bonnefoy                              *
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
 * @file ctfindvar.cpp
 * @brief Time variability search tool implementation
 * @author Simon Bonnefoy
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ctfindvar.hpp"
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include <cmath>

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#endif

/* __ Method name definitions ____________________________________________ */
#define G_INX2GTI                                  "ctfindvar::inx2gti(int&)"
#define G_FILL_CUBE                  "ctfindvar::fill_cube(GCTAObservation*)"
#define G_FILL_ALPHA_VECTOR       "ctfindvar::fill_alpha_vector(const int&, "\
                                                           "vector<double>&)"
#define G_INIT_CUBE                              "ctfindvar::init_cube(void)"
#define G_INIT_GTIS                              "ctfindvar::init_gtis(void)"

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
 *
 * Constructs empty time variability search tool.
 ***************************************************************************/
ctfindvar::ctfindvar(void) : ctobservation(CTFINDVAR_NAME, VERSION)
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
 * Constructs time variability search tool from an observation container.
 ***************************************************************************/
ctfindvar::ctfindvar(const GObservations& obs) :
           ctobservation(CTFINDVAR_NAME, VERSION, obs)
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
 * Constructs time variability search tool using command line arguments for
 * user parameter setting.
 ***************************************************************************/
ctfindvar::ctfindvar(int argc, char *argv[]) :
           ctobservation(CTFINDVAR_NAME, VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app search time variability tool.
 *
 * Constructs time variability search  tool from another time variability
 * search tool.
 ***************************************************************************/
ctfindvar::ctfindvar(const ctfindvar& app) : ctobservation(app)
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
 *
 * Destructs search time variability tool.
 ***************************************************************************/
ctfindvar::~ctfindvar(void)
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
 * @param[in] app Time variability search tool.
 * @return Time variability search tool.
 *
 * Assigns time variability search tool.
 ***************************************************************************/
ctfindvar& ctfindvar::operator=(const ctfindvar& app)
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
 * @brief Clear time variability search tool
 *
 * Clears time variability search tool.
 ***************************************************************************/
void ctfindvar::clear(void)
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
 * @brief Run time variability search tool
 ***************************************************************************/
void ctfindvar::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Create GTIs
    init_gtis();

    // Create counts cube
    create_cube();

    // Analyse cube
    analyse_cube();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save peak significance map and significance distributions
 *
 * Saves the peak significance and source signficance distributions.
 * Optionally, the generated counts cube can also be saved if 'outcube' has
 * been specified.
 ***************************************************************************/
void ctfindvar::save(void)
{
    // Write header
    log_header1(TERSE, "Saving results");

    // Write counts cube
    if ((*this)["outcube"].is_valid()) {

        // Get counts cube output filename
        GFilename outcube = (*this)["outcube"].filename();

        // Save counts cube
        m_counts.save(outcube, (*this)["clobber"].boolean());

        // Log saving
        log_value(TERSE, "Saving counts cube", outcube);
    }

    // Create a FITS file for storing the output
    GFits fits;

    // Write the most significant values for each pixel
    m_peaksigmap.write(fits, "PEAKSIGMAP");

    // Write source histograms
    write_source_histograms(fits);

    // Get output map filename
    GFilename outmap = (*this)["outmap"].filename();

    // Save output map
    fits.saveto(outmap, (*this)["clobber"].boolean());

    // Log saving
    log_value(TERSE, "Saving output map", outmap);

    // Save the model definition file
    if (m_model_above_thr.size() > 2) {

        // Get model output filename
        GFilename outmodel = (*this)["outmodel"].filename();

        // Save model
        m_model_above_thr.save(outmodel);

        // Log saving
        log_value(TERSE, "Saving output model", outmodel);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get the map index associated with a given time
 *
 * @param[in] time Time
 * @return Index of map in counts cube
 ***************************************************************************/
int ctfindvar::time2inx(const GTime& time)
{
    // Set default return value
    int map_index = -1;

    // Loop over all GTIs
    for (int i = 0; i < m_gti.size(); ++i) {

        // If interval contains the time then break
        if (m_gti[i].contains(time)) {
            map_index = i;
            break;
        }

    }

    // Return map index
    return map_index;
}


/***********************************************************************//**
 * @brief Return Good Time Interval for index
 *
 * @param[in] index Index
 * @return Good Time Interval
 ***************************************************************************/
GGti ctfindvar::inx2gti(const int& index)
{
    // Throw an exception if index is invalid
    if ((index < 0) || (index >= m_gti.size())) {
        throw GException::out_of_range(G_INX2GTI, "Time index", index, m_gti.size());
    }

    // Return Good Time Interval
    return m_gti[index];
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void ctfindvar::init_members(void)
{
    // Initialise members
    m_counts.clear();
    m_gti.clear();
    m_inmodel.clear();
    m_max_sig_dir.clear();
    m_minoff        = 0.0;
    m_sig_threshold = 0.0;
    m_peaksigmap.clear();
    m_pixsigsrc.clear();
    m_tstart.clear();
    m_tstop.clear();
    m_emin.clear();
    m_emax.clear();
    m_model_above_thr.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app search time variability tool.
 ***************************************************************************/
void ctfindvar::copy_members(const ctfindvar& app)
{
    // Copy members
    m_counts          = app.m_counts;
    m_gti             = app.m_gti;
    m_inmodel         = app.m_inmodel;
    m_max_sig_dir     = app.m_max_sig_dir;
    m_minoff          = app.m_minoff;
    m_sig_threshold   = app.m_sig_threshold;
    m_peaksigmap      = app.m_peaksigmap;
    m_pixsigsrc       = app.m_pixsigsrc;
    m_tstart          = app.m_tstart;
    m_tstop           = app.m_tstop;
    m_emin            = app.m_emin;
    m_emax            = app.m_emax;
    m_model_above_thr = app.m_model_above_thr;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctfindvar::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 ***************************************************************************/
void ctfindvar::get_parameters(void)
{
    // Load the observations
    setup_observations(m_obs, true, true, false);

    // Either load model or query source position parameters
    if ((*this)["inmodel"].is_valid()) {
        m_inmodel = GModels((*this)["inmodel"].filename());
    }
    else {
        (*this)["coordsys"].string();
        (*this)["xsrc"].real();
        (*this)["ysrc"].real();
    }

    // Query time limit parameters
    (*this)["tinterval"].real();
    if ((*this)["tmin"].is_valid()) {
        (*this)["tmin"].time();
    }
    if ((*this)["tmax"].is_valid()) {
        (*this)["tmax"].time();
    }

    // Create map to query spatial parameters
    create_map(m_obs);

    // Read energy limits
    m_emin = GEnergy((*this)["emin"].real(), "TeV");
    m_emax = GEnergy((*this)["emax"].real(), "TeV");

    // Get minimum counts for a bin to be considered in Noff calculation
    m_minoff = (*this)["minoff"].real();

    // Get minimum significance to set a source as variable
    m_sig_threshold = (*this)["threshold"].real();

    // Query smoothing parameters
    if ((*this)["smooth_kernel"].is_valid()) {
        (*this)["smooth_kernel"].string();
        (*this)["smooth_rad"].real();
    }

    // If needed later, query output filenames now
    if (read_ahead()) {
        (*this)["outcube"].query();
        (*this)["outmap"].query();
        (*this)["outmodel"].query();
    }

    // Set number of OpenMP threads
    #ifdef _OPENMP
    int nthreads = (*this)["nthreads"].integer();
    if (nthreads > 0) {
        omp_set_num_threads(nthreads);
    }
    #endif

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialize Good Time Intervals
 *
 * @exception GException::invalid_value
 *            Invalid start or stop time or insufficient time bins
 *
 * Initialise Good Time Intervals for the variability search.
 ***************************************************************************/
void ctfindvar::init_gtis(void)
{
    // Log header
    log_header1(NORMAL, "Initialise Good Time Intervals");

    // Clear time intervals
    m_gti.clear();

    // Set start and stop times
    m_tstart = get_tstart();
    m_tstop  = get_tstop();

    // Query time interval
    double tinterval = (*this)["tinterval"].real();

    // Compute time bins
    double tstart_sec = m_tstart.secs();
    double tstop_sec  = m_tstop.secs();
    int    bins       = (tstop_sec - tstart_sec) / tinterval + 0.5;

    // Log information
    log_value(NORMAL, "Reference time (mjd)", m_tstart.reference().mjdref());
    log_value(NORMAL, "Start time (sec)",    tstart_sec);
    log_value(NORMAL, "Stop time (sec)",     tstop_sec);
    log_value(NORMAL, "Total time (sec)",    tstop_sec-tstart_sec);
    log_value(NORMAL, "Time interval (sec)", tinterval);
    log_value(NORMAL, "Number of time bins", bins);

    // Throw exception if tstart == tstop
    if (m_tstart == m_tstop) {
        throw GException::invalid_value(G_INIT_GTIS,
                                        "Start time is equal to stop time");
    }
    else if (bins < 2) {
        std::string msg = "Method requires at least two time bins (" +
                          gammalib::str(bins) + " bins found).";
        throw GException::invalid_value(G_INIT_GTIS, msg);
    }

    // Set up Good Time Intervals by appending all intervals that overlap
    // with at least one of the observations
    for (int i = 0; i < bins; ++i) {

        // Get time interval
        GTime tstart(m_tstart + i*tinterval);
        GTime tstop(m_tstart + (i+1.0)*tinterval);
        GGti  gti(tstart, tstop);

        // Make sure there's an observation that overlaps with this GTI
        for (int k = 0; k < m_obs.size(); ++k) {

            // Get CTA observation
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[k]);

            // Continue only if observation is valid
            if (obs != NULL) {

                // Add GTI is the time intervals overlap
                if (gti_overlap(obs->gti(), gti) > 0.0) {
                    m_gti.push_back(gti);
                    break;
                }

            } // endif: observation was valid

        } // endfor: looped over observations

    } // endfor: looped over Good Time Intervals

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create counts cube
 *
 * Creates a counts cube comprising a sky map for each Good Time Interval.
 ***************************************************************************/
void ctfindvar::create_cube(void)
{
    // Log header
    log_header1(NORMAL, "Create counts cube");

    // Create the basic skymap
    m_counts = create_map(m_obs);

    // Resize to the appropriate number of time intervals
    m_counts.nmaps(m_gti.size());

    // Create the peaksigmap
    m_peaksigmap = m_counts.extract(1);

    // Loop over all unbinned CTA observations in the container
    #pragma omp parallel for
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get pointer to observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Fill the cube
        fill_cube(obs);

        // Dispose events to free memory
        obs->dispose_events();

    } // endfor: looped over observations

    // Smooth the maps if requested
    if ((*this)["smooth_kernel"].is_valid() &&
        (*this)["smooth_kernel"].string() != "NONE") {
        m_counts.smooth((*this)["smooth_kernel"].string(),
                        (*this)["smooth_rad"].real());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill the cube from the events in an observation
 *
 * @param[in] obs Pointer to CTA observation.
 *
 * Fills the cube from the events in an observation.
 ***************************************************************************/
void ctfindvar::fill_cube(GCTAObservation* obs)
{
    // Make sure that the observation holds a CTA event list. If this
    // is not the case then throw an exception.
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>
                                  (obs->events());
    if (events == NULL) {
        std::string msg = "CTA Observation does not contain an event "
                          "list. An event list is needed to fill the "
                          "counts cube.";
        throw GException::invalid_value(G_FILL_CUBE, msg);
    }

    // Get the RoI
    const GCTARoi& roi = events->roi();

    // Check for RoI sanity
    if (!roi.is_valid()) {
        std::string msg = "No valid RoI found in input observation "
                          "\""+obs->name()+"\". Run ctselect to specify "
                          "an RoI for this observation before running "
                          "ctbin.";
        throw GException::invalid_value(G_FILL_CUBE, msg);
    }

    // Initialise binning statistics
    int num_outside_roi  = 0;
    int num_invalid_wcs  = 0;
    int num_outside_ecut = 0;
    int num_outside_map  = 0;
    int num_outside_time = 0;
    int num_in_map       = 0;

    // Fill counts sky map
    for (int i = 0; i < events->size(); ++i) {

        // Get event
        const GCTAEventAtom* event = (*events)[i];

        // Check that the energy is valid
        if ((event->energy() < m_emin) || (event->energy() > m_emax)) {
            num_outside_ecut++;
            continue;
        }

        // Determine event time index
        int time_indx = time2inx(event->time());

        // Check if event is in the GTIs
        if ((time_indx > -1) && (time_indx < m_counts.nmaps())) {

            // Get event direction
            GSkyDir evnt_dir = GCTAInstDir(event->dir()).dir();

            // Skip event if it is outside the RoI
            if (roi.centre().dir().dist_deg(evnt_dir) > roi.radius()) {
                num_outside_roi++;
                continue;
            }

            // Determine sky pixel
            GSkyPixel pixel;
            try {
                pixel = m_counts.dir2pix(evnt_dir);
            }
            catch (std::exception &e) {
                num_invalid_wcs++;
                continue;
            }

            // Skip event if corresponding counts cube pixel is outside the
            // counts cube map range
            if (pixel.x() < -0.5 || pixel.x() > (m_counts.nx()-0.5) ||
                pixel.y() < -0.5 || pixel.y() > (m_counts.ny()-0.5)) {
                num_outside_map++;
                continue;
            }

            // Fill event in skymap
            #pragma omp critical(ctfindvar_fill_cube)
            m_counts(pixel, time_indx) += 1.0;

            // Increment number of maps
            num_in_map++;

        }

        // ... otherwise, event falls outside the specified time range
        else {
            num_outside_time++;
        }

    } // endfor: looped over all events

    // Log filling results
    #pragma omp critical(ctfindvar_fill_cube)
    {
        log_header3(TERSE, get_obs_header(obs));
        log_value(NORMAL, "Events in list", obs->events()->size());
        log_value(NORMAL, "Events in cube", num_in_map);
        log_value(NORMAL, "Events outside RoI", num_outside_roi);
        log_value(NORMAL, "Events outside energy", num_outside_ecut);
        log_value(NORMAL, "Events with invalid WCS", num_invalid_wcs);
        log_value(NORMAL, "Events outside cube area", num_outside_map);
        log_value(NORMAL, "Events outside time bins", num_outside_time);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Analyse all pixels of counts cube
 ***************************************************************************/
void ctfindvar::analyse_cube(void)
{
    // Log header
    log_header1(NORMAL, "Analyse sky map pixels");

    // Get number of GTIs
    const int nbins = m_gti.size();

    // Get relevant sky map pixels
    std::vector<int> srcInxPix = get_pixels();

    // Preparing the histogram with significance evolution for the source
    // center and the highest sig-pix
    m_pixsigsrc = GNdarray(srcInxPix.size(), nbins);
    m_pixsigmax = GNdarray(nbins);

    // Initialise maximum significance
    double max_sig = 0;

    // looping over all the pixels in the cube
    #pragma omp parallel for
    for (int ipix = 0; ipix < m_counts.npix(); ++ipix) {

        // Compute total counts for pixel
        double total_counts = 0.0;
        for (int k = 0; k < nbins; ++k) {
            total_counts += m_counts(ipix, k);
        }

        // Getting the variability significance for the current pixel
        GNdarray pixSig = get_variability_sig(ipix);

        // Getting the significance evolution for the source
        for (int src = 0; src < srcInxPix.size(); ++src) {

            // Store the distribution if the source is located at this position
            if (srcInxPix[src] == ipix) {

                #pragma omp critical(ctfindvar_analyse_cube)
                for (int i = 0; i < nbins; i++) {
                    m_pixsigsrc(src,i) = pixSig(i);
                }
            }

        } // endfor: looped over sources

        // Make sure that this code is not executed in parallel
        #pragma omp critical(ctfindvar_analyse_cube)
        {  
            // Determine the pixel with the highest significance
            if (max(pixSig) > max_sig) {
                max_sig       = max(pixSig);
                m_max_sig_dir = m_counts.inx2dir(ipix);
                m_pixsigmax   = pixSig;
            }

            // Storing sig in skymap
            m_peaksigmap(ipix) = max(pixSig);

            // If the significance is greater than the significance threshold,
            // the position is stored in a model
            if (max(pixSig) > m_sig_threshold) {
                GSkyDir dir_pix = m_counts.inx2dir(ipix);
                m_model_above_thr.append(sky_model(dir_pix));
            }
        }

    } // endfor: looped over pixels

    // Log results
    log_value(NORMAL, "Maximum significance", max_sig);
    log_value(NORMAL, "Right Ascension of maximum", m_max_sig_dir.ra_deg());
    log_value(NORMAL, "Declination of maximum", m_max_sig_dir.dec_deg());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return pixel vector
 *
 * @return Vector of sky map pixels
 *
 * Returns a vector of sky map pixels that corresponds to the source model
 * positions.
 ***************************************************************************/
std::vector<int> ctfindvar::get_pixels(void)
{
    // Initialise pixel array
    std::vector<int> srcInxPix;

    // If there are input models then extract
    if (m_inmodel.size() > 0) {

        // Loop over models
        for (int i = 0; i < m_inmodel.size(); ++i) {

            // Extract the model
            GModelSky* model = dynamic_cast<GModelSky*>(m_inmodel[i]);

            // If model is a sky model then extract the model centre. Note
            // that this only works for a model that results a GSkyRegionCircle
            // object as region.
            if (model != NULL) {

                // Get the source position
                GModelSpatial*    spatial = model->spatial();
                GSkyRegion*       region  = spatial->region();
                GSkyRegionCircle* circle  = dynamic_cast<GSkyRegionCircle*>(region);

                // If region circle is valid then
                if (circle != NULL) {
                    int ipix = m_counts.dir2inx(circle->centre());
                    srcInxPix.push_back(ipix);
                    continue;
                }

            } // endif: model was a sky model

            // If we are still alive then remove the model since it is neither
            // a sky model, nor does it have a region circle
            m_inmodel.remove(i);
            i--;

        } // endfor: looped over all models

    } // endif: there were input models

    // ... otherwise query the source position
    else {

        // Query source direction
        GSkyDir dir;
        if ((*this)["coordsys"].string() == "CEL") {
            dir.radec_deg((*this)["xsrc"].real(),(*this)["ysrc"].real());
        }
        else {
            dir.lb_deg((*this)["xsrc"].real(),(*this)["ysrc"].real());
        }

        // Store the index of this source
        int ipix = m_counts.dir2inx(dir);
        srcInxPix.push_back(ipix);

    } // endelse: queried source position

    // Return pixel array
    return srcInxPix;
}


/***********************************************************************//**
 * @brief Return sky model for a given sky direction
 *
 * @param[in]  dir Sky direction
 * @return Sky model for a given sky direction
 *
 * Returns a point source sky model at the position of the current pixel
 * with a power law spectral model corresponding to 1% of the Crab flux.
 *
 * The model name is built from the celestial coordinates of the source.
 ***************************************************************************/
GModelSky ctfindvar::sky_model(const GSkyDir& dir) const
{
    // Set point source spatial component from sky direction
    GModelSpatialPointSource spatial(dir);

    // Set fake spectral component (powerlaw with 1% of Crab)
    GModelSpectralPlaw spectral(5.7e-18, -2.48, GEnergy(0.3, "TeV"));

    // Create a sky model
    GModelSky model(spatial, spectral);

    // Build model name from celestial coordinates
    double      ra   = dir.ra_deg();
    double      dec  = dir.dec_deg();
    std::string name = gammalib::str(ra, 4);
    if (dec < 0.0) {
        name += "-";
    }
    else {
        name += "+";
    }
    name += gammalib::str(std::abs(dec), 4);

    // Set model name
    model.name(name);

    // Return model
    return model;
}


/***********************************************************************//**
 * @brief Get the significance of variability for a given skymap pixel
 * 
 * @param[in] ipix Sky map pixel index.
 * @return Histogram of significances.
 * 
 * Significance is computed according to Li & Ma equation 17:
 *    see: https://doi.org/10.1086/161295
 * Some modification is made in order to handle the case where Non = 0, 
 * but Noff != 0.
 ***************************************************************************/
GNdarray ctfindvar::get_variability_sig(const int& ipix)
{
    // Initialise result
    GNdarray sig_histogram(m_gti.size());

    // Initialise vectors
    std::vector<bool> accepted_bin_bckg_vector(m_gti.size(), true);

    // Get alpha values for specified pixel
    std::vector<double> alphas = get_alphas(ipix);

    // Exclude all bins from the background estimate for which the
    // number of events is below "minoff" or for which alpha is zero
    for (int i = 0; i < m_gti.size(); ++i) {
        if (m_counts(ipix, i) < m_minoff || alphas[i] == 0.0) {
            accepted_bin_bckg_vector[i] = false;
        }
    }

    // Loop over pixels until all background pixels are removed
    bool background_validated = false;
    while (background_validated == false) {

        // Signal that background was validated
        background_validated = true;

        // Loop over all the GTIs of the pixel
        for (int i = 0; i < m_gti.size(); ++i) {

            // The GTI is discarded from background calculation and not
            // checked again
            if (!accepted_bin_bckg_vector[i]) {
                continue;
            }

            // ... otherwise GTI is selected, and we loop over all the
            // others
            double noff  = 0.0;
            double alpha = 0.0;
            for (int j = 0; j < m_gti.size(); ++j) {
                if (j != i && accepted_bin_bckg_vector[j]) {
                    noff  += m_counts(ipix, j);
                    alpha += alphas[j];
                }
            }

            // The background is averaged on the number of bins -1
            double non = m_counts(ipix, i);
            alpha      = alphas[i]/alpha;

            // Compute sensitivity in Gaussian sigma using Li & Ma
            double alpha1 = alpha + 1.0;
            double ntotal = non+noff;
            double arg1   = non/ntotal;
            double arg2   = noff/ntotal;
            double sig    = 0.0;
            if (alpha != 0.0) {
                if (noff == 0.0) {
                    sig = 0.0;
                }
                else if (non > 0.0) {
                    double term1  = non  * std::log((alpha1/alpha)*arg1);
                    double term2  = noff * std::log(alpha1*arg2);
                    sig  = std::sqrt(2.0 * (term1 + term2));
                }
                else {
                    sig = std::sqrt(2.0 * noff * std::log(alpha1));
                }
                sig *= (non < alpha*noff) ? -1.0 : 1.0;
            }

            // Update the significance
            sig_histogram(i) = sig;

            // If bin is significant, it's removed from the background and
            // we loop again
            if (sig > m_sig_threshold) {
                accepted_bin_bckg_vector[i] = false;
                background_validated        = false;
            }

        } // endfor: looped over all GTIs

    } // endwhile: iterate until converged

    // Return significance histogram
    return sig_histogram;
}


/***********************************************************************//**
 * @brief Get alpha vector
 *
 * @param[in] ipix Sky map pixel index
 * @return Vector containing effective exposure
 ***************************************************************************/
std::vector<double> ctfindvar::get_alphas(const int& ipix) const
{
    // Initialise alpha vector
    std::vector<double> alphas(m_gti.size(), 0.0);

    // Convert pixel index into sky direction
    GSkyDir pix_dir = m_counts.inx2dir(ipix);

    // Loop over all observations
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        const GCTAObservation* obs = dynamic_cast<const GCTAObservation*>(m_obs[i]);

        // Fall through if observation is not valid
        if (obs == NULL) {
            continue;
        }

        // Get RoI. Fall through if  observation does not overlap with this
        // pixel position
        GCTARoi roi        = obs->roi();
        GSkyDir roi_centre = roi.centre().dir();
        if (roi_centre.dist_deg(pix_dir) > roi.radius()) {
            continue;
        }

        // Convert sky direction to instrument direction
        GCTAInstDir instdir = obs->pointing().instdir(pix_dir);

        // Extract the good time intervals of the observation
        const GGti& obs_gti(obs->gti());

        // Loop over all time bins
        for (int j = 0; j < m_gti.size(); ++j) {

            // Make sure observation overlaps with this time interval
            double exposure = gti_overlap(m_gti[j], obs_gti);
            if (exposure > 0.0) {

                // Get IRF response
                const GCTAResponseIrf* rsp =
                      dynamic_cast<const GCTAResponseIrf*>(obs->response());

                // Throw an exception if observation has no instrument response function
                if (rsp == NULL) {
                    std::string msg = 
                        "No response information available for "+
                        get_obs_header(obs)+" to compute IRF background. "
                        "Please specify response information or use "
                        "another background subtraction method.";
                    throw GException::invalid_value(G_FILL_ALPHA_VECTOR, msg);
                }

                // Get IRF background template
                const GCTABackground* bkg = rsp->background();

                // Throw an exception if observation has no IRF background template
                if (bkg == NULL) {
                    std::string msg = 
                        "No IRF background template found in instrument "
                        "response function for "+
                        get_obs_header(obs)+". Please specify an instrument "
                        "response function containing a background template.";
                    throw GException::invalid_value(G_FILL_ALPHA_VECTOR, msg);
                }

                // Add up background rate
                #pragma omp critical(ctfindvar_get_alphas) 
                alphas[j] += exposure * bkg->rate_ebin(instdir, m_emin, m_emax);

            } // endif: exposure was positive

        } // endfor: looped over all time bins

    } // endfor: looped over all observations

    // Return alphas
    return alphas;
} 


/***********************************************************************//**
 * @brief Returns number of seconds that two GTIs overlap
 *
 * @param[in] gti1      First GTI
 * @param[in] gti2      Second GTI
 * @return Number of seconds that @p gti1 and @p gti2 overlap
 ***************************************************************************/
double ctfindvar::gti_overlap(const GGti& gti1, const GGti& gti2) const
{
    // Initialise no overlap
    double overlap = 0.0;

    // Set first and second GTI so that gti_1st starts before gti_2nd
    const GGti& gti_1st = (gti1.tstart() <= gti2.tstart()) ? gti1 : gti2;
    const GGti& gti_2nd = (gti1.tstart()  > gti2.tstart()) ? gti1 : gti2;

    // Compute overlap
    if (gti_1st.tstop() > gti_2nd.tstart()) {
        if (gti_1st.tstop() <= gti_2nd.tstop()) {
            overlap = gti_1st.tstop() - gti_2nd.tstart();
        }
        else if (gti_1st.tstop() > gti_2nd.tstop()) {
            overlap = gti_2nd.tstop() - gti_2nd.tstart();
        }
    }

    // Return overlap
    return overlap;
}


/***********************************************************************//**
 * @brief Get start time
 *
 * @return Start time
 ***************************************************************************/
GTime ctfindvar::get_tstart(void)
{
    // Initialise start time with large unreasonable Julian date
    GTime tstart("JD 999999999");

    // If tmin parameter is valid then set the start time
    if ((*this)["tmin"].is_valid()) {

        // Query tmin parameter
        tstart = (*this)["tmin"].time();

        // Log action
        log_string(NORMAL, "Setting start time \"tmin\" parameter");

    } // endif: tmin parameter was valid

    // ... otherwise set start time from observations
    else {

        // Loop over all observations
        for (int k = 0; k < m_obs.size(); ++k) {

            // Get CTA observation
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[k]);

            // Continue only if observation is valid
            if (obs != NULL) {

                // Update start time if the start time of the observation is
                // smaller than the current start time
                if (obs->gti().tstart() < tstart) {
                    tstart = obs->gti().tstart();
                }

            } // endif: observation was valid

        } // endfor: looped over observations

        // Log action
        log_string(NORMAL, "Setting start time from observations");

    } // endelse: no valid time specified

    // Return start time
    return tstart;
}


/***********************************************************************//**
 * @brief Get stop time
 *
 * @return Stop time
 ***************************************************************************/
GTime ctfindvar::get_tstop(void)
{
    // Initialise stop time with small unreasonable Julian date
    GTime tstop("JD 0");

    // If tmax parameter is valid then set the stop time
    if ((*this)["tmax"].is_valid()) {

        // Query tmax parameter
        tstop = (*this)["tmax"].time();

        // Log action
        log_string(NORMAL, "Setting stop time \"tmax\" parameter");

    } // endif: tmax parameter was valid

    // ... otherwise set stop time from observations
    else {

        // Loop over all observations
        for (int k = 0; k < m_obs.size(); ++k) {

            // Get CTA observation
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[k]);

            // Continue only if observation is valid
            if (obs != NULL) {

                // Update stop time if the stop time of the observation is
                // larger than the current stop time
                if (obs->gti().tstop() > tstop) {
                    tstop = obs->gti().tstop();
                }

            } // endif: observation was valid

        } // endfor: looped over observations

        // Log action
        log_string(NORMAL, "Setting stop time from observations");

    } // endelse: no valid time specified

    // Return stop time
    return tstop;
}


/***********************************************************************//**
 * @brief Write source histograms in FITS format
 * 
 * @param[in,out] fits FITS file
 *
 * This method creates two FITS binary table extensions and appends them to
 * the FITS file.
 *
 * The first extension is named "SOURCE SIGNIFICANCE" and contains the
 * significance of the variability as function of time bin for the pixel
 * with the maximum significance and for each source.
 *
 * The second extension is named "SOURCE POSITION" and contains the Right
 * Ascension and Declination for the maximum significance pixel and for each
 * source.
 ***************************************************************************/
void ctfindvar::write_source_histograms(GFits& fits)
{
    // Setup data storage for time and source significances
    double   tinterval = (*this)["tinterval"].real();
    int      nbins     = (m_tstop-m_tstart) / tinterval + 0.5;
    int      nsrc      = m_pixsigsrc.shape()[0];
    GNdarray time_info(2, nbins);
    GNdarray max_pixel_info(nbins);
    GNdarray src_info(nsrc, nbins);

    // Loop over time bins
    for (int i = 0; i < nbins; ++i) {

        // Set time bin mid point
        GTime tstart_bin(m_tstart + i * tinterval);
        GTime tstop_bin(m_tstart + (i+1) * tinterval);
        GTime midpoint((tstop_bin.secs()+tstart_bin.secs())*0.5);

        // Store the time of this bin
        time_info(0,i) = tstart_bin.mjd();
        time_info(1,i) = tstop_bin.mjd();

        // Get the index associated with this time
        int index = time2inx(midpoint);

        // Loop over every source
        for (int k = 0; k < nsrc; ++k) {
            src_info(k,i) = (index >= 0) ? m_pixsigsrc(k,index) : 0.0;
        }

        // Fill maximum pixel information
        max_pixel_info(i) = (index >= 0) ? m_pixsigmax(index) : 0.0;

    } // endfor: looped over time bins

    // Create binary tables
    GFitsBinTable table_signif(nbins);
    GFitsBinTable table_pos(nsrc+1);
    table_signif.extname("SOURCE SIGNIFICANCE");
    table_pos.extname("SOURCE POSITION");

    // Create columns for time bins
    GFitsTableDoubleCol time_start("TSTART", nbins);
    GFitsTableDoubleCol time_stop("TSTOP", nbins);
    GFitsTableDoubleCol time_sigma("MAXSIGPIXEL", nbins);
    time_start.unit("MJD");
    time_stop.unit("MJD");
    time_sigma.unit("sigma");

    // Create columns for source position
    GFitsTableStringCol pos_name("SOURCE", nsrc+1, 20);
    GFitsTableDoubleCol pos_ra("RA", nsrc+1);
    GFitsTableDoubleCol pos_dec("DEC", nsrc+1);
    pos_ra.unit("deg");
    pos_dec.unit("deg");

    // Fill time bin columns
    for (int i = 0; i < nbins; ++i) {
        time_start(i) = time_info(0,i);
        time_stop(i)  = time_info(1,i);
        time_sigma(i) = max_pixel_info(i);
    }

    // Fill source columns
    pos_name(0) = "MAXSIGPIXEL";
    pos_ra(0)   = m_max_sig_dir.ra_deg();
    pos_dec(0)  = m_max_sig_dir.dec_deg();

    // Append columns to tables
    table_signif.append(time_start);
    table_signif.append(time_stop);
    table_signif.append(time_sigma);

    // Append the distribution for each source
    for (int k = 0; k < nsrc; ++k) {

        // Initialise source name and direction
        std::string src_name("SOURCE");
        GSkyDir     src_dir;

        // Get the name from the input model file
        if (m_inmodel.size() > k) {

            // Extract the model
            GModelSky* model = dynamic_cast<GModelSky*>(m_inmodel[k]);

            // If model is a sky model then extract the model centre. Note
            // that this only works for a model that results a GSkyRegionCircle
            // object as region.
            if (model != NULL) {

                // Get the source position
                GModelSpatial*    spatial = model->spatial();
                GSkyRegion*       region  = spatial->region();
                GSkyRegionCircle* circle  = dynamic_cast<GSkyRegionCircle*>(region);

                // If region circle is valid then extract source name and
                // sky direction
                if (circle != NULL) {
                    src_name = m_inmodel[k]->name();
                    src_dir  = circle->centre();
                }

            } // endif: model was a sky model

        } // endif: there were enough models in the input model

        // ... otherwise set the position from the xsrc and ysrc parameters
        else {
            if ((*this)["coordsys"].string() == "CEL") {
                src_dir.radec_deg((*this)["xsrc"].real(), (*this)["ysrc"].real());
            }
            else {
                src_dir.lb_deg((*this)["xsrc"].real(), (*this)["ysrc"].real());
            }
        }

        // Create new time column for source
        GFitsTableDoubleCol column(src_name, nbins);
        column.unit("sigma");

        // Fill column
        for (int i = 0; i < nbins; ++i) {
            column(i) = src_info(k, i);
        }

        // Append column to time table
        table_signif.append(column);

        // Store source infirmation
        pos_name(k+1) = src_name;
        pos_ra(k+1)   = src_dir.ra_deg();
        pos_dec(k+1)  = src_dir.dec_deg();

    } // endfor: looped over sources

    // Append columns to tables
    table_pos.append(pos_name);
    table_pos.append(pos_ra);
    table_pos.append(pos_dec);

    // Append tables to FITS file
    fits.append(table_signif);
    fits.append(table_pos);

    // Return
    return;
}
