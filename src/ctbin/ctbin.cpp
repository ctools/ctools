/***************************************************************************
 *                        ctbin - Event binning tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @file ctbin.cpp
 * @brief Event binning tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctbin.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_FILL_CUBE                      "ctbin::fill_cube(GCTAObservation*)"
#define G_GET_PARAMETERS                            "ctbin::get_parameters()"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_PLAW_WEIGHTS      //!< Use power law in energy weight computation

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty ctbin tool.
 ***************************************************************************/
ctbin::ctbin(void) : ctool(CTBIN_NAME, CTBIN_VERSION)
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
 * Constructs ctbin tool from an observation container.
 ***************************************************************************/
ctbin::ctbin(const GObservations& obs) : ctool(CTBIN_NAME, CTBIN_VERSION)
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
 *
 * Constructs ctbin tool using command line arguments for user parameter
 * setting.
 ***************************************************************************/
ctbin::ctbin(int argc, char *argv[]) : 
       ctool(CTBIN_NAME, CTBIN_VERSION, argc, argv)
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
 *
 * Constructs ctbin tool from another ctbin instance.
 ***************************************************************************/
ctbin::ctbin(const ctbin& app) : ctool(app)
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
 * Destructs ctbin tool.
 ***************************************************************************/
ctbin::~ctbin(void)
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
 * @param[in] app ctbin tool.
 * @return ctbin tool.
 *
 * Assigns ctbin tool.
 ***************************************************************************/
ctbin& ctbin::operator=(const ctbin& app)
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
 * @brief Clear ctbin tool
 *
 * Clears ctbin tool.
 ***************************************************************************/
void ctbin::clear(void)
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
 * @brief Run the ctbin tool
 *
 * Gets the user parameters and loops over all CTA observations in the
 * observation container to bin the events into a single counts cube. All
 * observations in the observation container that do not contain CTA event
 * lists will be skipped.
 ***************************************************************************/
void ctbin::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        log.header1(gammalib::number("Observation", m_obs.size()));
        log << m_obs << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1(gammalib::number("Bin observation", m_obs.size()));
    }

    // Loop over all observations in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Write header for observation
        if (logTerse()) {
            log.header3(get_obs_header(m_obs[i]));
        }

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Skip observation if it's not CTA
        if (obs == NULL) {
            if (logTerse()) {
                log << " Skipping ";
                log << m_obs[i]->instrument();
                log << " observation" << std::endl;
            }
            continue;
        }

        // Skip observation if we have a binned observation
        if (obs->eventtype() == "CountsCube") {
            if (logTerse()) {
                log << " Skipping binned ";
                log << obs->instrument();
                log << " observation" << std::endl;
            }
            continue;
        }

        // Fill the cube
        fill_cube(obs);

        // Dispose events to free memory
        obs->dispose_events();

    } // endfor: looped over observations

    // Normalize weights counts map
    if (m_livetime > 0.0) {
        for (int iebin = 0; iebin < m_weights.nmaps(); ++iebin) {
            for (int index = 0; index < m_weights.npix(); ++index) {
                m_weights(index, iebin) /= m_livetime;
            }
        }
    }

    // Set a single cube in the observation container
    obs_cube();

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        log.header1(gammalib::number("Binned observation", m_obs.size()));
        log << m_obs << std::endl;
    }

    // Optionally publish counts cube
    if (m_publish) {
        publish();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save counts cube
 *
 * Saves the counts cube into a FITS file.
 ***************************************************************************/
void ctbin::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save counts cube");
    }

    // Get counts cube filename
    m_outcube = (*this)["outcube"].filename();

    // Save only if filename is non-empty
    if (!m_outcube.is_empty()) {

        // Get CTA observation from observation container
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

        // Save only if observation is valid
        if (obs != NULL) {
        
            // Log filename
            if (logTerse()) {
                log << gammalib::parformat("Counts cube file");
                log << m_outcube.url() << std::endl;
            }
            
            // Save cube
            obs->save(m_outcube, clobber());

        } // endif: observation was valid

    } // endif: outcube file was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Publish counts cube
 *
 * @param[in] name Counts cube name.
 ***************************************************************************/
void ctbin::publish(const std::string& name)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Publish counts cube");
    }

    // Set default name is user name is empty
    std::string user_name(name);
    if (user_name.empty()) {
        user_name = CTBIN_NAME;
    }

    // Log filename
    if (logTerse()) {
        log << gammalib::parformat("Counts cube name");
        log << user_name << std::endl;
    }

    // Publish counts cube
    m_counts.publish(user_name);

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
void ctbin::init_members(void)
{
    // Initialise members
    m_outcube.clear();
    m_usepnt  = false;
    m_publish = false;

    // Initialise protected members
    m_obs.clear();
    m_counts.clear();
    m_weights.clear();
    m_ebounds.clear();
    m_gti.clear();
    m_ontime   = 0.0;
    m_livetime = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctbin::copy_members(const ctbin& app)
{
    // Copy attributes
    m_outcube = app.m_outcube;
    m_usepnt  = app.m_usepnt;
    m_publish = app.m_publish;

    // Copy protected members
    m_obs      = app.m_obs;
    m_counts   = app.m_counts;
    m_weights  = app.m_weights;
    m_ebounds  = app.m_ebounds;
    m_gti      = app.m_gti;
    m_ontime   = app.m_ontime;
    m_livetime = app.m_livetime;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctbin::free_members(void)
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
void ctbin::get_parameters(void)
{
    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {

        // Throw exception if no input observation file is given
        require_inobs(G_GET_PARAMETERS);

        // Throw exception if counts cube is given
        require_inobs_nocube(G_GET_PARAMETERS);

        // Get observation container without response (not needed)
        m_obs = get_observations(false);

    } // endif: there was no observation in the container

    // Create an event cube based on task parameters
    GCTAEventCube cube = create_cube(m_obs);

    // Get the skymap from the cube and initialise all counts cube bins and
    // weights to zero
    m_counts  = cube.counts();
    m_counts  = 0.0;
    m_weights = m_counts;

    // Get energy boundaries
    m_ebounds  = cube.ebounds();

    // Get remaining parameters
    m_publish = (*this)["publish"].boolean();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outcube = (*this)["outcube"].filename();
    }

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill events into counts cube
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No event list or valid RoI found in observation.
 *
 * Fills the events from an event list in the counts cube setup by init_cube.
 ***************************************************************************/
void ctbin::fill_cube(GCTAObservation* obs)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Make sure that the observation holds a CTA event list. If this
        // is not the case then throw an exception.
        const GCTAEventList* events =
              dynamic_cast<const GCTAEventList*>(obs->events());
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
        int num_outside_map  = 0;
        int num_outside_ebds = 0;
        int num_in_map       = 0;

        // Fill counts sky map
        for (int i = 0; i < events->size(); ++i) {

            // Get event
            const GCTAEventAtom* event = (*events)[i];

            // Determine sky pixel
            GCTAInstDir* inst  = (GCTAInstDir*)&(event->dir());
            GSkyDir      dir   = inst->dir();
            GSkyPixel    pixel = m_counts.dir2pix(dir);

            // Skip if pixel is outside RoI
            if (roi.centre().dir().dist_deg(dir) > roi.radius()) {
                num_outside_roi++;
                continue;
            }

            // Skip if pixel is out of map range
            if (pixel.x() < -0.5 || pixel.x() > (m_counts.nx()-0.5) ||
                pixel.y() < -0.5 || pixel.y() > (m_counts.ny()-0.5)) {
                num_outside_map++;
                continue;
            }

            // Determine energy bin
            int iebin = m_ebounds.index(event->energy());

            // Skip if energy bin is outside counts cube boundaries
            if (iebin == -1) {
                num_outside_ebds++;
                continue;
            }

            // Fill event in skymap
            m_counts(pixel, iebin) += 1.0;
            num_in_map++;

        } // endfor: looped over all events

        // Get RoI and energy weights
        GSkyMap             roi_weights    = set_roi_weights(roi);
        std::vector<double> energy_weights = set_energy_weights(events->ebounds());

        // Add weights multiplied by livetime to weight map. The map is later
        // divided by the cumulative livetime in run().
        for (int iebin = 0; iebin < energy_weights.size(); ++iebin) {
            for (int index = 0; index < roi_weights.npix(); ++index) {
                m_weights(index, iebin) += roi_weights(index) *
                                           energy_weights[iebin] *
                                           obs->livetime();
            }
        }

        // Append GTIs
        m_gti.extend(events->gti());

        // Update ontime and livetime
        m_ontime   += obs->ontime();
        m_livetime += obs->livetime();

        // Log filling results
        if (logTerse()) {
            log << gammalib::parformat("Events in list");
            log << obs->events()->size() << std::endl;
            log << gammalib::parformat("Events in cube");
            log << num_in_map << std::endl;
            log << gammalib::parformat("Event bins outside RoI");
            log << num_outside_roi << std::endl;
            log << gammalib::parformat("Events outside cube area");
            log << num_outside_map << std::endl;
            log << gammalib::parformat("Events outside energy bins");
            log << num_outside_ebds << std::endl;
        }

        // Log cube
        if (logExplicit()) {
            log.header1("Counts cube");
            log << m_counts << std::endl;
        }

    } // endif: observation was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create output observation container.
 *
 * Creates an output observation container that combines all input CTA
 * observation into a single stacked observation. All non-CTA observations
 * and all binned CTA observations that were present in the observation
 * container are append to the observation container so that they can
 * be used by other tools. The method furthermore conserves any response
 * information in case that a single CTA observation is provided to support
 * binned analysis.
 ***************************************************************************/
void ctbin::obs_cube(void)
{
    // If we have only a single CTA observation in the container, then
    // keep that observation and just attach the event cube to it. Reset
    // the filename, otherwise we still will have the old event filename
    // in the log file.
    if (m_obs.size() == 1) {

        // Attach event cube to CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);
        if (obs != NULL) {

            // Change the event type if we had an unbinned observation
            if (obs->eventtype() == "EventList") {

                // Assign cube to the observation
                obs->events(this->cube());

            }

            // ... otherwise the input observation was binned and hence
            // skipped. In that case we simply append an empty counts
            // cube
            else {

                // Create empty counts cube
                GCTAEventCube cube(m_counts, m_weights, m_ebounds, obs->gti());

                // Assign empty cube to observation
                obs->events(cube);

            }

            // Reset file name
            obs->eventfile("");

        } // endif: obervation was valid

    } // endif: we only had one observation in the container

    // ... otherwise put a single CTA observation in container and
    // append all observations that have not been used in the binning
    // to the container
    else {

        // Allocate observation container
        GObservations container;

        // Allocate CTA observation.
        GCTAObservation obs;

        // Attach event cube to CTA observation
        obs.events(this->cube());

        // Compute average pointing direction for all CTA event lists
        double ra     = 0.0;
        double dec    = 0.0;
        double number = 0.0;
        for (int i = 0; i < m_obs.size(); ++i) {
            GCTAObservation* cta = dynamic_cast<GCTAObservation*>(m_obs[i]);
            if ((cta != NULL) && (cta->eventtype() == "EventList")) {
                ra     += cta->pointing().dir().ra();
                dec    += cta->pointing().dir().dec();
                number += 1.0;
            }
        }
        if (number > 0.0) {
            ra  /= number;
            dec /= number;
        }
        GSkyDir dir;
        dir.radec(ra, dec);
        GCTAPointing pointing(dir);

        // Compute deadtime correction
        double deadc = (m_ontime > 0.0) ? m_livetime / m_ontime : 0.0;

        // Set CTA observation attributes
        obs.pointing(pointing);
        obs.obs_id(0);
        obs.ra_obj(dir.ra_deg());   //!< Dummy
        obs.dec_obj(dir.dec_deg()); //!< Dummy
        obs.ontime(m_ontime);
        obs.livetime(m_livetime);
        obs.deadc(deadc);

        // Set models in observation container
        container.models(m_obs.models());

        // Append CTA observation
        container.append(obs);
        
        // Copy over all remaining non-CTA observations
        for (int i = 0; i < m_obs.size(); ++i) {
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
            if (obs == NULL) {
                container.append(*m_obs[i]);
            }
            else if (obs->eventtype() != "EventList") {
                container.append(*m_obs[i]);
            }
        }

        // Set observation container
        m_obs = container;

    } // endelse: there was not a single CTA observation

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set RoI weights.
 *
 * @param[in] roi Region of Interest.
 * @return Sky map with Region of Interest weights.
 *
 * Computes a sky map of Region of Interest (RoI) weights. All pixels of the
 * sky map that are fully comprised within the RoI will be set to unity.
 * Pixels that are only partly comprised will be devided into a 10x10 subgrid,
 * and a weight corresponding to the fraction of subgrid points that are
 * fully comprised within the RoI will be assigned to the sky map pixel.
 * The precision of the weights is about 1%.
 ***************************************************************************/
GSkyMap ctbin::set_roi_weights(const GCTARoi& roi) const
{
    // Set number of points in subgrid
    const int n_subgrid = 10;

    // Initialise ROI weight sky map with zero weights
    GSkyMap weights(m_counts.extract(0));
    weights = 0.0;

    // Compute save RoI radius. All pixels with bin centre within that radius
    // will be set to unity. In case that the counts cube is not a WCS map
    // no pixels will be considered within the safe radius.
    double      margin(roi.radius());
    const GWcs* wcs(dynamic_cast<const GWcs*>(m_counts.projection()));
    if (wcs != NULL) {
        double dx = wcs->cdelt(0);
        double dy = wcs->cdelt(1);
        margin    = (dx > dy) ? dx : dy;
    }
    double safe_radius = roi.radius() - margin;

    // Loop over all sky map pixels
    for (int i = 0; i < weights.npix(); ++i) {

        // Get pixel sky direction
        GSkyDir dir(weights.inx2dir(i));

        // If pixel is within the safe radius then set the pixel weight to
        // unity ...
        if (roi.centre().dir().dist_deg(dir) <= safe_radius) {
            weights(i) = 1.0;
        }

        // ... otherwise define a subgrid a check all positions within that
        // subgrid that fall within the RoI. This is an approximation that
        // determines the weight of the pixels
        else {

            // Set pixel centre
            GSkyPixel centre(weights.inx2pix(i));

            // Compute subgrid pixel step size
            double step = 1.0 / double(n_subgrid-1);

            // Compute sky direction weight
            double weight = 1.0 / (double(n_subgrid) * double(n_subgrid));

            // Loop over X subgrid
            for (int ix = 0; ix < n_subgrid; ++ix) {

                // Compute X pixel value
                double x = ix * step - 0.5;

                // Loop over Y subgrid
                for (int iy = 0; iy < n_subgrid; ++iy) {

                    // Compute Y pixel value
                    double y = iy * step - 0.5;

                    // Set pixel value
                    GSkyPixel pixel(centre.x()+x, centre.y()+y);
                    
                    // Get pixel sky direction
                    GSkyDir dir(weights.pix2dir(pixel));

                    // If sky direction is within the RoI then add the direction
                    // to the weight for this pixel
                    if (roi.centre().dir().dist_deg(dir) <= roi.radius()) {
                        weights(i) += weight;
                    }

                } // endfor: loop over Y subgrid

            } // endfor: loop over X subgrid

        } // endelse: pixel was outside safe radius

    } // endfor: looped over all pixel sky maps

    // Return
    return weights;
}


/***********************************************************************//**
 * @brief Set energy weights.
 *
 * @param[in] ebounds Energy boundaries.
 * @return Vector of energy boundary weights.
 *
 * Computes a vector of energy boundary weights from the overlap of the
 * event list energy boundaries with the counts cube boundaries. It is
 * assumed that within a bin the events are distributed according to a power
 * law with spectral index of -2, and the weights are computed according
 * to the fraction of this power law that is comprised within a counts cube
 * bin.
 *
 * This method properly handles event lists with multiple energy boundaries.
 ***************************************************************************/
std::vector<double> ctbin::set_energy_weights(const GEbounds& ebounds) const
{
    // Initialise energy boundary weights with zero weights
    std::vector<double> weights(m_ebounds.size(), 0.0);

    // Loop over all energy bins of counts cube
    for (int i = 0; i < m_ebounds.size(); ++i) {

        // Get energy boundaries of counts cube bin
        GEnergy emin(m_ebounds.emin(i));
        GEnergy emax(m_ebounds.emax(i));

        // Loop over all energy boundaries of the event list and add the
        // coverage fraction to the weights
        for (int k = 0; k < ebounds.size(); ++k) {

            // If there is no overlap then check next boundary
            if ((ebounds.emax(k) <= emin) || (ebounds.emin(k) >= emax)) {
                continue;
            }

            // Get minimum and maximum energy of overlap region
            GEnergy elow  = (ebounds.emin(k) > emin) ? ebounds.emin(k) : emin;
            GEnergy ehigh = (ebounds.emax(k) < emax) ? ebounds.emax(k) : emax;

            // Compute overlap weight
            #if defined(G_PLAW_WEIGHTS)
            double w_overlap = gammalib::plaw_photon_flux(elow.TeV(),
                                                          ehigh.TeV(),
                                                          1.0, -2.0);
            double w_full    = gammalib::plaw_photon_flux(emin.TeV(),
                                                          emax.TeV(),
                                                          1.0, -2.0);
            double weight = w_overlap / w_full;
            #else
            double weight = (ehigh - elow) / (emax - emin);
            #endif

            // Add overlap weight
            weights[i] += weight;

        } // endfor: looped over all energy boundaries of the event list

    } // endfor: looped over all energy boundaries of the counts cube

    // Return
    return weights;
}
