/***************************************************************************
 *                ctobssim - CTA observation simulator tool                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file ctobssim.cpp
 * @brief CTA observation simulator tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctobssim.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_SIMULATE_SOURCE       "ctobssim::simulate_source(GCTAObservation*,"\
                                                                " GPhotons&)"

/* __ Constants __________________________________________________________ */
const double g_roi_margin = 0.5;      //!< Simulation radius margin (degrees)

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
ctobssim::ctobssim(void) : GApplication(CTOBSSIM_NAME, CTOBSSIM_VERSION)
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observations constructor
 *
 * @param[in] obs Observation container.
 *
 * Constructs application by initialising from an observation container.
 ***************************************************************************/
ctobssim::ctobssim(GObservations obs) :
          GApplication(CTOBSSIM_NAME, CTOBSSIM_VERSION)
{
    // Initialise members
    init_members();

    // Set observations
    m_obs = obs;

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Command line constructor
 *
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 *
 * Constructs application using command line arguments for parameter setting.
 ***************************************************************************/
ctobssim::ctobssim(int argc, char *argv[]) : 
          GApplication(CTOBSSIM_NAME, CTOBSSIM_VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Application.
 ***************************************************************************/
ctobssim::ctobssim(const ctobssim& app) : GApplication(app)
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
ctobssim::~ctobssim(void)
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
ctobssim& ctobssim::operator= (const ctobssim& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->GApplication::operator=(app);

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
 * @brief Clear application
 ***************************************************************************/
void ctobssim::clear(void)
{
    // Free members
    free_members();
    this->GApplication::free_members();

    // Initialise members
    this->GApplication::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Execute application
 *
 * This method runs the simulation and saves the results.
 ***************************************************************************/
void ctobssim::execute(void)
{
    // Signal that some parameters should be read ahead
    m_read_ahead = true;

    // Run the simulation
    run();

    // Save results
    save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Simulate event data
 *
 * This method runs the simulation. Results are not saved by this method.
 * Invoke "save" to save the results.
 ***************************************************************************/
void ctobssim::run(void)
{
    // Switch screen logging on in debug mode
    if (logDebug()) {
        log.cout(true);
    }

    // Get parameters
    get_parameters();

    // Write input parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Write seed values into logger
    if (logTerse()) {
        log << std::endl;
        log.header1("Seed values");
        for (int i = 0; i < m_rans.size(); ++i) {
            log << gammalib::parformat("Seed "+gammalib::str(i));
            log << gammalib::str(m_rans[i].seed()) << std::endl;
        }
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

    // Write header
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Simulate observations");
        }
        else {
            log.header1("Simulate observation");
        }
    }

    // Initialise counters
    int n_observations = 0;

    // From here on the code can be parallelized if OpenMP support
    // is enabled. The code in the following block corresponds to the
    // code that will be executed in each thread
    #pragma omp parallel
    {
        // Each thread will have it's own logger to avoid conflicts
        GLog wrklog;
        if (logDebug()) {
            wrklog.cout(true);
        }

        // Copy configuration from application logger to thread logger
        wrklog.date(log.date());
        wrklog.name(log.name());

        // Set a big value to avoid flushing
        wrklog.max_size(10000000);

        // Loop over all observation in the container. If OpenMP support
        // is enabled, this looped will be parallelized.
        #pragma omp for
        for (int i = 0; i < m_obs.size(); ++i) {

            // Get CTA observation
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

            // Continue only if observation is a CTA observation
            if (obs != NULL) {

                // Write header for observation
                if (logTerse()) {
                    if (obs->name().length() > 1) {
                        wrklog.header3("Observation "+obs->name());
                    }
                    else {
                        wrklog.header3("Observation");
                    }
                }

                // Increment counter
                n_observations++;

                // Simulate source events
                simulate_source(obs, m_obs.models(), m_rans[i], &wrklog);

                // Simulate source events
                simulate_background(obs, m_obs.models(), m_rans[i], &wrklog);

            } // endif: CTA observation found

        } // endfor: looped over observations

        // At the end, the content of the thread logger is added to
        // the application logger
        #pragma omp critical (log)
        {
            log << wrklog;
        }

    } // end pragma omp parallel
    
    // If more than a single observation has been handled then make sure that
    // an XML file will be used for storage
    if (n_observations > 1) {
        m_use_xml = true;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save the selected event list(s)
 *
 * This method saves the selected event list(s) into FITS file(s). There are
 * two modes, depending on the m_use_xml flag.
 *
 * If m_use_xml is true, all selected event list(s) will be saved into FITS
 * files, where the output filenames are constructued from the input
 * filenames by prepending the m_prefix string to name. Any path information
 * will be stripped form the input name, hence event files will be written
 * into the local working directory (unless some path information is present
 * in the prefix). In addition, an XML file will be created that gathers
 * the filename information for the selected event list(s). If an XML file
 * was present on input, all metadata information will be copied from this
 * input file.
 *
 * If m_use_xml is false, the selected event list will be saved into a FITS
 * file.
 ***************************************************************************/
void ctobssim::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Save observations");
        }
        else {
            log.header1("Save observation");
        }
    }

    // Case A: Save event file(s) and XML metadata information
    if (m_use_xml) {
        save_xml();
    }

    // Case B: Save event file as FITS file
    else {
        save_fits();
    }

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
void ctobssim::get_parameters(void)
{
    // If there are no observations in container then add a single CTA
    // observation using the parameters from the parameter file
    if (m_obs.size() == 0) {

        // Get CTA observation parameters
        m_infile = (*this)["infile"].filename();
        m_caldb  = (*this)["caldb"].string();
        m_irf    = (*this)["irf"].string();
        m_ra     = (*this)["ra"].real();
        m_dec    = (*this)["dec"].real();

        // Set pointing direction
        GCTAPointing pnt;
        GSkyDir      skydir;
        skydir.radec_deg(m_ra, m_dec);
        pnt.dir(skydir);

        // Allocate CTA observation
        GCTAObservation obs;

        // Set CTA observation attributes
        obs.pointing(pnt);
        obs.response(m_irf, m_caldb);

        // Set event list (queries remaining parameters)
        set_list(&obs);

        // Append CTA observation to container
        m_obs.append(obs);

        // Load models into container
        m_obs.models(m_infile);

        // Signal that no XML file should be used for storage
        m_use_xml = false;

    } // endif: there was no observation in the container

    // ... otherwise we signal that an XML file should be written
    else {
        m_use_xml = true;
    }

    // Get other parameters
    m_seed = (*this)["seed"].integer();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (m_read_ahead) {
        m_outfile = (*this)["outfile"].filename();
        m_prefix  = (*this)["prefix"].string();
    }

    // Initialise random number generators. We initialise here one random
    // number generator per observation so that each observation will
    // get it's own random number generator. This will lead to identical
    // results independently of code parallelization with OpenMP. The
    // seeds for all random number generators are derived randomly but
    // fully deterministacally from the seed parameter, so that a given
    // seed parameter leads always to the same set of simulated events, and
    // this independently of parallelization.

    // Get a random number generator for seed determination
    GRan master(m_seed);

    // Allocate vector of random number generator seeds
    std::vector<unsigned long long int> seeds;

    // Loop over all observations in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Allocate new seed value
        unsigned long long int new_seed;

        // Determine new seed value. We make sure that the new seed
        // value has not been used already for another observation.
        bool repeat = false;
        do {
            new_seed = (unsigned long long int)(master.uniform() * 1.0e10) +
                       m_seed;
            repeat   = false;
            for (int j = 0; j < seeds.size(); ++j) {
                if (new_seed == seeds[j]) {
                    repeat = true;
                    break;
                }
            }
        } while(repeat);

        // Add the seed to the vector for bookkeeping
        seeds.push_back(new_seed);

        // Use the seed to create a random number generator for the
        // actual observation
        m_rans.push_back(GRan(new_seed));

    } // endfor: looped over observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set empty CTA event list
 *
 * @param[in] obs CTA observation.
 *
 * Attaches an empty event list to CTA observation. The method also sets the
 * pointing direction using the m_ra and m_dec members, the ROI based on
 * m_ra, m_dec and m_rad, a single GTI based on m_tmin and m_tmax, and a
 * single energy boundary based on m_emin and m_emax. The method furthermore
 * sets the ontime, livetime and deadtime correction factor.
 ***************************************************************************/
void ctobssim::set_list(GCTAObservation* obs)
{
    // Continue only if observation is valid
    if (obs != NULL) {

        // Get CTA observation parameters
        m_ra    = (*this)["ra"].real();
        m_dec   = (*this)["dec"].real();
        m_rad   = (*this)["rad"].real();
        m_tmin  = (*this)["tmin"].real();
        m_tmax  = (*this)["tmax"].real();
        m_emin  = (*this)["emin"].real();
        m_emax  = (*this)["emax"].real();
        m_deadc = (*this)["deadc"].real();

        // Allocate CTA event list
        GCTAEventList events;

        // Set pointing direction
        GCTAPointing pnt;
        GSkyDir      skydir;
        skydir.radec_deg(m_ra, m_dec);
        pnt.dir(skydir);

        // Set ROI
        GCTAInstDir instdir(skydir);
        GCTARoi     roi(instdir, m_rad);

        // Set GTI
        GGti  gti(m_cta_ref);
        GTime tstart;
        GTime tstop;
        tstart.set(m_tmin, m_cta_ref);
        tstop.set(m_tmax, m_cta_ref);
        gti.append(tstart, tstop);

        // Set energy boundaries
        GEbounds ebounds;
        GEnergy  emin;
        GEnergy  emax;
        emin.TeV(m_emin);
        emax.TeV(m_emax);
        ebounds.append(emin, emax);

        // Set CTA event list attributes
        events.roi(roi);
        events.gti(gti);
        events.ebounds(ebounds);

        // Attach event list to CTA observation
        obs->events(events);

        // Set observation ontime, livetime and deadtime correction factor
        obs->ontime(gti.ontime());
        obs->livetime(gti.ontime()*m_deadc);
        obs->deadc(m_deadc);

    } // endif: oberservation was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Simulate source events from photon list
 *
 * @param[in] obs Pointer on CTA observation.
 * @param[in] models Model list.
 * @param[in] ran Random number generator.
 * @param[in] wrklog Pointer to logger.
 *
 * Simulate source events from a photon list for a given CTA observation.
 * The events are stored in as event list in the observation.
 *
 * This method does nothing if the observation pointer is NULL. It also
 * verifies if the observation has a valid pointing and response.
 ***************************************************************************/
void ctobssim::simulate_source(GCTAObservation* obs, const GModels& models,
                               GRan& ran, GLog* wrklog)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // If no logger is specified then use the default logger
        if (wrklog == NULL) {
            wrklog = &log;
        }

        // Get pointer on CTA response. Throw an exception if the response
        // is not defined.
        const GCTAResponse& rsp = obs->response();

        // Make sure that the observation holds a CTA event list. If this
        // is not the case then allocate and attach a CTA event list now.
        if (dynamic_cast<const GCTAEventList*>(obs->events()) == NULL) {
            set_list(obs);
        }

        // Get pointer on event list (circumvent const correctness)
        GCTAEventList* events = static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));

        // Extract simulation region.
        GSkyDir dir = events->roi().centre().dir();
        double  rad = events->roi().radius() + g_roi_margin;

        // Dump simulation cone information
        if (logNormal()) {
            *wrklog << gammalib::parformat("Simulation area");
            *wrklog << m_area << " cm2" << std::endl;
            *wrklog << gammalib::parformat("Simulation cone");
            *wrklog << "RA=" << dir.ra_deg() << " deg";
            *wrklog << ", Dec=" << dir.dec_deg() << " deg";
            *wrklog << ", r=" << rad << " deg" << std::endl;
        }

        // Initialise indentation for logging
        int indent = 0;

        // Loop over all Good Time Intervals
        for (int it = 0; it < events->gti().size(); ++it) {

            // Extract time interval
            GTime tmin = events->gti().tstart(it);
            GTime tmax = events->gti().tstop(it);

            // Dump time interval
            if (logNormal()) {
                if (events->gti().size() > 1) {
                    indent++;
                    wrklog->indent(indent);
                }
                *wrklog << gammalib::parformat("Time interval", indent);
                *wrklog << tmin.convert(m_cta_ref);
                *wrklog << " - ";
                *wrklog << tmax.convert(m_cta_ref);
                *wrklog << " s" << std::endl;
            }

            // Loop over all energy boundaries
            for (int ie = 0; ie <  events->ebounds().size(); ++ie) {

                // Extract energy boundaries
                GEnergy emin = events->ebounds().emin(ie);
                GEnergy emax = events->ebounds().emax(ie);

                // Dump energy range
                if (logNormal()) {
                    if (events->ebounds().size() > 1) {
                        indent++;
                        wrklog->indent(indent);
                    }
                    *wrklog << gammalib::parformat("Energy range", indent);
                    *wrklog << emin << " - " << emax << std::endl;
                }

                // Loop over all sky models
                for (int i = 0; i < models.size(); ++i) {

                    // Get sky model (NULL if not a sky model)
                    const GModelSky* model =
                          dynamic_cast<const GModelSky*>(models[i]);

                    // If we have a sky model then simulate photons
                    if (model != NULL) {

                        // Determine duration of a time slice by limiting the
                        // number of simulated photons to m_max_photons.
                        // The photon rate is estimated from the model flux
                        // and used to set the duration of the time slice.
                        double flux     = model->spectral()->flux(emin, emax);
                        double rate     = flux * m_area;
                        double duration = 1800.0;  // default: 1800 sec
                        if (rate > 0.0) {
                            duration = m_max_photons / rate;
                            if (duration < 1.0) {  // not <1 sec
                                duration = 1.0;
                            }
                            else if (duration > 180000.0) { // not >50 hr
                                duration = 180000.0;
                            }
                        }
                        GTime tslice(duration, "sec");

                        // To reduce memory requirements we split long time
                        // intervals into several slices.
                        GTime tstart = tmin;
                        GTime tstop  = tstart + tslice;

                        // Initialise cumulative photon counters
                        int nphotons = 0;

                        // Loop over time slices
                        while (tstart < tmax) {

                            // Make sure that tstop <= tmax
                            if (tstop > tmax) {
                                tstop = tmax;
                            }

                            // Dump time slice
                            if (logExplicit()) {
                                if (tmax - tmin > tslice) {
                                    indent++;
                                    wrklog->indent(indent);
                                }
                                *wrklog << gammalib::parformat("Time slice", indent);
                                *wrklog << tstart.convert(m_cta_ref);
                                *wrklog << " - ";
                                *wrklog << tstop.convert(m_cta_ref);
                                *wrklog << " s" << std::endl;
                            }

                            // Get photons
                            GPhotons photons = model->mc(m_area, dir, rad,
                                                         emin, emax,
                                                         tstart, tstop, ran);

                            // Dump number of simulated photons
                            if (logExplicit()) {
                                *wrklog << gammalib::parformat("MC source photons/slice", indent);
                                *wrklog << photons.size();
                                if (model->name().length() > 0) {
                                    *wrklog << " [" << model->name() << "]";
                                }
                                *wrklog << std::endl;
                            }

                            // Simulate events from photons
                            for (int i = 0; i < photons.size(); ++i) {

                                // Increment photon counter
                                nphotons++;

                                // Simulate event. Note that this method
                                // includes the deadtime correction.
                                GCTAEventAtom* event = rsp.mc(m_area,
                                                              photons[i],
                                                              *obs,
                                                              ran);
                                if (event != NULL) {

                                    // Use event only if it falls within ROI
                                    if (events->roi().contains(*event)) {
                                        event->event_id(m_event_id);
                                        events->append(*event);
                                        m_event_id++;
                                    }
                                    delete event;
                                }

                            } // endfor: looped over events

                            // Go to next time slice
                            tstart = tstop;
                            tstop  = tstart + tslice;

                            // Reset indentation
                            if (logExplicit()) {
                                if (tmax - tmin > tslice) {
                                    indent--;
                                    wrklog->indent(indent);
                                }
                            }

                        } // endwhile: looped over time slices

                        // Dump simulation results
                        if (logNormal()) {
                            *wrklog << gammalib::parformat("MC source photons", indent);
                            *wrklog << nphotons;
                            if (model->name().length() > 0) {
                                *wrklog << " [" << model->name() << "]";
                            }
                            *wrklog << std::endl;
                            *wrklog << gammalib::parformat("MC source events", indent);
                            *wrklog << events->size();
                            if (model->name().length() > 0) {
                                *wrklog << " [" << model->name() << "]";
                            }
                            *wrklog << std::endl;

                        }

                    } // endif: model was a sky model

                } // endfor: looped over models

                // Dump simulation results
                if (logNormal()) {
                    *wrklog << gammalib::parformat("MC source events", indent);
                    *wrklog << events->size();
                    *wrklog << " (all source models)";
                    *wrklog << std::endl;
                }

                // Reset indentation
                if (logNormal()) {
                    if (events->ebounds().size() > 1) {
                        indent--;
                        wrklog->indent(indent);

                    }
                }

            } // endfor: looped over all energy boundaries

            // Reset indentation
            if (logNormal()) {
                if (events->gti().size() > 1) {
                    indent--;
                    wrklog->indent(indent);
                }
            }

        } // endfor: looped over all time intervals

        // Reset indentation
        wrklog->indent(0);

    } // endif: observation pointer was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Simulate background events from model
 *
 * @param[in] obs Pointer on CTA observation.
 * @param[in] models Models.
 * @param[in] ran Random number generator.
 * @param[in] wrklog Pointer to logger.
 *
 * Simulate background events from models. The events are stored as event
 * list in the observation.
 *
 * This method does nothing if the observation pointer is NULL.
 ***************************************************************************/
void ctobssim::simulate_background(GCTAObservation* obs,
                                   const GModels&   models,
                                   GRan&            ran, 
                                   GLog*            wrklog)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // If no logger is specified then use the default logger
        if (wrklog == NULL) {
            wrklog = &log;
        }

        // Make sure that the observation holds a CTA event list. If this
        // is not the case then allocate and attach a CTA event list now.
        if (dynamic_cast<const GCTAEventList*>(obs->events()) == NULL) {
            set_list(obs);
        }

        // Get pointer on event list (circumvent const correctness)
        GCTAEventList* events =
            static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));

        // Loop over all models
        for (int i = 0; i < models.size(); ++i) {

            // Get data model (NULL if not a data model)
            const GModelData* model =
                dynamic_cast<const GModelData*>(models[i]);

            // If we have a data model that applies to CTA then simulate
            // events
            if (model != NULL && model->isvalid("CTA", "")) {

                // Get simulated CTA event list. Note that this method
                // includes the deadtime correction.
                GCTAEventList* list =
                    dynamic_cast<GCTAEventList*>(model->mc(*obs, ran));

                // Continue only if we got a CTA event list
                if (list != NULL) {

                    // Reserves space for events
                    events->reserve(list->size()+events->size());

                    // Append events
                    for (int k = 0; k < list->size(); k++) {

                        // Get event pointer
                        GCTAEventAtom* event = (*list)[k];

                        // Use event only if it falls within ROI
                        if (events->roi().contains(*event)) {

                            // Set event identifier
                            event->event_id(m_event_id);
                            m_event_id++;

                            // Append event
                            events->append(*event);

                        } // endif: event was within ROI

                    } // endfor: looped over all events

                    // Dump simulation results
                    if (logNormal()) {
                        *wrklog << gammalib::parformat("MC background events");
                        *wrklog << list->size() << std::endl;
                    }

                    // Free event list
                    delete list;

                } // endif: we had a CTA event list

            } // endif: model was valid

        } // endfor: looped over all models

    } // endif: observation pointer was valid

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
 *
 * The thrown area is fixed to pi*(2500^2) m2, which is the same value that
 * is used in the Monte Carlo simulations (information from Konrad Bernloehr)
 ***************************************************************************/
void ctobssim::init_members(void)
{
    // Initialise user parameters
    m_infile.clear();
    m_outfile.clear();
    m_prefix.clear();
    m_caldb.clear();
    m_irf.clear();
    m_seed  =   1;
    m_ra    = 0.0;
    m_dec   = 0.0;
    m_rad   = 0.0;
    m_tmin  = 0.0;
    m_tmax  = 0.0;
    m_emin  = 0.0;
    m_emax  = 0.0;
    m_deadc = 1.0;

    // Initialise protected members
    m_rans.clear();
    m_obs.clear();
    m_use_xml    = false;
    m_read_ahead = false;

    // Set fixed parameters
    m_area        = 19634954.0 * 1.0e4; //!< pi*(2500^2) m^2
    m_max_photons = 1000000;            //!< Maximum number of photons / time slice

    // Set CTA time reference. G_CTA_MJDREF is the CTA reference MJD,
    // which is defined in GCTALib.hpp. This is somehow a kluge. We need
    // a better mechanism to implement the CTA reference MJD.
    m_cta_ref.set(G_CTA_MJDREF, "s", "TT", "LOCAL");

    // Initialise first event identifier
    m_event_id = 1;

    // Set logger properties
    log.date(true);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctobssim::copy_members(const ctobssim& app)
{
    // Copy user parameters
    m_infile   = app.m_infile;
    m_outfile  = app.m_outfile;
    m_prefix   = app.m_prefix;
    m_caldb    = app.m_caldb;
    m_irf      = app.m_irf;
    m_seed     = app.m_seed;
    m_ra       = app.m_ra;
    m_dec      = app.m_dec;
    m_rad      = app.m_rad;
    m_tmin     = app.m_tmin;
    m_tmax     = app.m_tmax;
    m_emin     = app.m_emin;
    m_emax     = app.m_emax;
    m_deadc    = app.m_deadc;

    // Copy protected members
    m_area        = app.m_area;
    m_max_photons = app.m_max_photons;
    m_rans        = app.m_rans;
    m_obs         = app.m_obs;
    m_use_xml     = app.m_use_xml;
    m_read_ahead  = app.m_read_ahead;
    m_cta_ref     = app.m_cta_ref;
    m_event_id    = app.m_event_id;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctobssim::free_members(void)
{
    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save event list in FITS format.
 *
 * Save the event list as a FITS file. The filename of the FITS file is
 * specified by the outfile parameter.
 ***************************************************************************/
void ctobssim::save_fits(void)
{
    // Get output filename
    m_outfile = (*this)["outfile"].filename();

    // Get CTA observation from observation container
    GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

    // Save observation into FITS file
    obs->save(m_outfile, clobber());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save event list(s) in XML format.
 *
 * Save the event list(s) into FITS files and write the file path information
 * into a XML file. The filename of the XML file is specified by the outfile
 * parameter, the filename(s) of the event lists are built by prepending a
 * prefix to the input event list filenames. Any path present in the input
 * filename will be stripped, i.e. the event list(s) will be written in the
 * local working directory (unless a path is specified in the prefix).
 ***************************************************************************/
void ctobssim::save_xml(void)
{
    // Get output filename and prefix
    m_outfile = (*this)["outfile"].filename();
    m_prefix  = (*this)["prefix"].string();

    // Initialise file number
    int file_num = 0;

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Handle only CTA observations
        if (obs != NULL) {

            // Set event output file name
            std::string outfile = m_prefix + gammalib::str(file_num) + ".fits";

            // Store output file name in observation
            obs->eventfile(outfile);

            // Save observation into FITS file
            obs->save(outfile, clobber());

            // Increment file number
            file_num++;


        } // endif: observation was a CTA observations

    } // endfor: looped over observations

    // Save observations in XML file
    m_obs.save(m_outfile);

    // Return
    return;
}
