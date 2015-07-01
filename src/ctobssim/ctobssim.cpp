/***************************************************************************
 *                  ctobssim - Observation simulator tool                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2015 by Juergen Knoedlseder                         *
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
 * @brief Observation simulator tool implementation
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
#define G_SIMULATE_SOURCE      "ctobssim::simulate_source(GCTAObservation*, "\
                                                    "GModels&, GRan&, GLog*)"

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
ctobssim::ctobssim(void) : ctool(CTOBSSIM_NAME, CTOBSSIM_VERSION)
{
    // Initialise members
    init_members();

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
ctobssim::ctobssim(const GObservations& obs) :
          ctool(CTOBSSIM_NAME, CTOBSSIM_VERSION)
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
 * Constructs application using command line arguments for parameter setting.
 ***************************************************************************/
ctobssim::ctobssim(int argc, char *argv[]) :
          ctool(CTOBSSIM_NAME, CTOBSSIM_VERSION, argc, argv)
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
ctobssim::ctobssim(const ctobssim& app) : ctool(app)
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
ctobssim& ctobssim::operator=(const ctobssim& app)
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
 * @brief Clear application
 ***************************************************************************/
void ctobssim::clear(void)
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

    // Special mode: if read ahead is specified we know that we called
    // the execute() method, hence files are saved immediately and event
    // lists are disposed afterwards.
    if (read_ahead()) {
        m_save_and_dispose = true;
    }

    // Determine the number of valid CTA observations, set energy dispersion flag
    // for all CTA observations and save old values in save_edisp vector
    int               n_observations = 0;
    std::vector<bool> save_edisp;
    save_edisp.assign(m_obs.size(), false);
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            save_edisp[i] = obs->response()->apply_edisp();
            obs->response()->apply_edisp(m_apply_edisp);
            n_observations++;
        }
    }

    // If more than a single observation has been handled then make sure that
    // an XML file will be used for storage
    if (n_observations > 1) {
        m_use_xml = true;
    }

    // Write execution mode into logger
    if (logTerse()) {
        log << std::endl;
        log.header1("Execution mode");
        log << gammalib::parformat("Event list management");
        if (m_save_and_dispose) {
            log << "Save and dispose (reduces memory needs)" << std::endl;
        }
        else {
            log << "Keep events in memory" << std::endl;
        }
        log << gammalib::parformat("Output format");
        if (m_use_xml) {
            log << "Write Observation Definition XML file" << std::endl;
        }
        else {
            log << "Write single event list FITS file" << std::endl;
        }
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

        // Allocate and initialize copies for multi-threading
        GModels models(m_obs.models());

        // Copy configuration from application logger to thread logger
        wrklog.date(log.date());
        wrklog.name(log.name());

        // Set a big value to avoid flushing
        wrklog.max_size(10000000);

        // Loop over all observation in the container. If OpenMP support
        // is enabled, this loop will be parallelized.
        #pragma omp for
        for (int i = 0; i < m_obs.size(); ++i) {

            // Get pointer on CTA observation
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

                // Work on a clone of the CTA observation. This makes sure that
                // any memory allocated for computing (for example a response
                // cache) is properly de-allocated on exit of this run
                GCTAObservation obs_clone = *obs;

                // Save number of events before entering simulation
                int events_before = obs_clone.events()->size();

                // Simulate source events
                simulate_source(&obs_clone, models, m_rans[i], &wrklog);

                // Simulate source events
                simulate_background(&obs_clone, models, m_rans[i], &wrklog);

                // Dump simulation results
                if (logNormal()) {
                    wrklog << gammalib::parformat("MC events");
                    wrklog << obs_clone.events()->size() - events_before;
                    wrklog << " (all models)";
                    wrklog << std::endl;
                }

                // Append the event list to the original observation
                obs->events(*(obs_clone.events()));

                // If requested, event lists are saved immediately
                if (m_save_and_dispose) {

                    // Set event output file name. If multiple observations are
                    // handled, build the filename from prefix and observation
                    // index. Otherwise use the outfile parameter.
                    std::string outfile;
                    if (m_use_xml) {
                        m_prefix = (*this)["prefix"].string();
                        outfile  = m_prefix + gammalib::str(i) + ".fits";
                    }
                    else {
                        outfile  = (*this)["outevents"].filename();
                    }

                    // Store output file name in original observation
                    obs->eventfile(outfile);

                    // Save observation into FITS file. This is a critical zone
                    // to avoid multiple threads writing simultaneously
                    #pragma omp critical
                    {
                        obs_clone.save(outfile, clobber());
                    }

                    // Dispose events
                    obs->dispose_events();

                }

                // ... otherwise append the event list to the original observation
                /*
                else {
                    obs->events(*(obs_clone.events()));
                }
                */

            } // endif: CTA observation found

        } // endfor: looped over observations

        // At the end, the content of the thread logger is added to
        // the application logger
        #pragma omp critical (log)
        {
            log << wrklog;
        }

    } // end pragma omp parallel

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
    m_outevents.clear();
    m_prefix.clear();
    m_seed        = 1;
    m_apply_edisp = false;
    m_max_rate    = 1.0e6;

    // Initialise protected members
    m_rans.clear();
    m_obs.clear();
    m_save_and_dispose = false;

    // Set fixed parameters
    m_area        = 19634954.0 * 1.0e4; //!< pi*(2500^2) m^2
    m_max_photons = 1000000;            //!< Maximum number of photons / time slice

    // Initialise first event identifier
    m_event_id = 1;

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
    m_outevents   = app.m_outevents;
    m_prefix      = app.m_prefix;
    m_seed        = app.m_seed;
    m_apply_edisp = app.m_apply_edisp;
    m_max_rate    = app.m_max_rate;

    // Copy protected members
    m_area             = app.m_area;
    m_max_photons      = app.m_max_photons;
    m_rans             = app.m_rans;
    m_obs              = app.m_obs;
    m_save_and_dispose = app.m_save_and_dispose;
    m_event_id         = app.m_event_id;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctobssim::free_members(void)
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
void ctobssim::get_parameters(void)
{
    // Initialise seed vector
    m_rans.clear();

    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {
        m_obs = get_observations();
    }

    // ... otherwise make sure that observation boundaries are set
    else {
        set_obs_bounds(m_obs);
    }

    // Read model definition file if required
    if (m_obs.models().size() == 0) {
       std::string inmodel = (*this)["inmodel"].filename();
       m_obs.models(inmodel);
    }

    // Get other parameters
    m_seed        = (*this)["seed"].integer();
    m_apply_edisp = (*this)["edisp"].boolean();
    m_max_rate    = (*this)["maxrate"].real();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outevents = (*this)["outevents"].filename();
        m_prefix    = (*this)["prefix"].string();
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

        // Get CTA response
        const GCTAResponseIrf* rsp =
            dynamic_cast<const GCTAResponseIrf*>(obs->response());
        if (rsp == NULL) {
            std::string msg = "Response is not an IRF response.\n" +
                              obs->response()->print();
            throw GException::invalid_value(G_SIMULATE_SOURCE, msg);
        }

        // Get pointer on event list (circumvent const correctness)
        GCTAEventList* events =
            static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));

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

                // Set true photon energy limits for simulation. If observation
                // has energy dispersion then add margin
                GEnergy e_true_min = emin;
                GEnergy e_true_max = emax;
                if (rsp->use_edisp()) {
                    e_true_min = rsp->ebounds(e_true_min).emin();
                    e_true_max = rsp->ebounds(e_true_max).emax();
                }

                // Dump energy range
                if (logNormal()) {
                    if (events->ebounds().size() > 1) {
                        indent++;
                        wrklog->indent(indent);
                    }
                    *wrklog << gammalib::parformat("Photon energy range", indent);
                    *wrklog << e_true_min << " - " << e_true_max << std::endl;
                    *wrklog << gammalib::parformat("Event energy range", indent);
                    *wrklog << emin << " - " << emax << std::endl;
                }

                // Loop over all sky models
                for (int i = 0; i < models.size(); ++i) {

                    // Get sky model (NULL if not a sky model)
                    const GModelSky* model =
                          dynamic_cast<const GModelSky*>(models[i]);

                    // If we have a sky model that applies to the present
                    // observation then simulate events
                    if (model != NULL &&
                        model->is_valid(obs->instrument(), obs->id())) {

                        // Determine duration of a time slice by limiting the
                        // number of simulated photons to m_max_photons.
                        // The photon rate is estimated from the model flux
                        // and used to set the duration of the time slice.
                        double flux     = get_model_flux(model, emin, emax,
                                                         dir, rad);
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

                        // Dump photon rate
                        if (logNormal()) {
                            *wrklog << gammalib::parformat("Photon rate", indent);
                            *wrklog << rate << " photons/sec";
                            if (model->name().length() > 0) {
                                *wrklog << " [" << model->name() << "]";
                            }
                            *wrklog << std::endl;
                        }

                        // If photon rate exceeds the maximum photon rate
                        // then throw an exception
                        if (rate > m_max_rate) {
                            std::string modnam = (model->name().length() > 0) ?
                                                 model->name() : "Unknown";
                            std::string msg    = "Photon rate "+
                                                 gammalib::str(rate)+
                                                 " photons/sec for model \""+
                                                 modnam+"\" exceeds maximum"
                                                 " allowed photon rate of "+
                                                 gammalib::str(m_max_rate)+
                                                 " photons/sec. Please check"
                                                 " the model parameters for"
                                                 " model \""+modnam+"\" or"
                                                 " increase the value of the"
                                                 " hidden \"maxrate\""
                                                 " parameter.";
                            throw GException::invalid_value(G_SIMULATE_SOURCE, msg);
                        }

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
                                *wrklog << " s";
                                if (model->name().length() > 0) {
                                    *wrklog << " [" << model->name() << "]";
                                }
                                *wrklog << std::endl;
                            }

                            // Get photons
                            GPhotons photons = model->mc(m_area, dir, rad,
                                                         e_true_min, e_true_max,
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
                                GCTAEventAtom* event = rsp->mc(m_area,
                                                               photons[i],
                                                               *obs,
                                                               ran);

                                if (event != NULL) {

                                    // Use event only if it falls within ROI
                                    // energy boundary and time slice
                                    if (events->roi().contains(*event) &&
                                        event->energy() >= emin &&
                                        event->energy() <= emax &&
                                        event->time() >= tstart &&
                                        event->time() <= tstop) {
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
 * @brief Determine sky model flux
 *
 * @param[in] model Sky model.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 * @param[in] centre Centre of region for photon rate determination.
 * @param[in] radius Radius of region for photon rate determination.
 * @return Model flux (photons/cm2/sec).
 ***************************************************************************/
double ctobssim::get_model_flux(const GModelSky* model,
                                const GEnergy&   emin,
                                const GEnergy&   emax,
                                const GSkyDir&   centre,
                                const double&    radius)
{
    // Initialise de-allocation flag
    bool free_spectral = false;

    // Get pointer to spectral model
    const GModelSpectral* spectral = model->spectral();

    // If the spatial model is a diffuse cube then create a node function
    // spectral model that is the product of the diffuse cube node function
    // and the spectral model evaluated at the energies of the node function
    GModelSpatialDiffuseCube* cube = dynamic_cast<GModelSpatialDiffuseCube*>(model->spatial());
    if (cube != NULL) {

        // Set MC cone
        cube->set_mc_cone(centre, radius);

        // Allocate node function to replace the spectral component
        GModelSpectralNodes* nodes = new GModelSpectralNodes(cube->spectrum());
        for (int i = 0; i < nodes->nodes(); ++i) {
            GEnergy energy    = nodes->energy(i);
            GTime   time;                              // Dummy time
            double  intensity = nodes->intensity(i);
            double  norm      = spectral->eval(energy, time);
            nodes->intensity(i, norm*intensity);
        }

        // Signal that node function needs to be de-allocated later
        free_spectral = true;

        // Set the spectral model pointer to the node function
        spectral = nodes;

    } // endif: spatial model was a diffuse cube
    
    // Compute flux within [emin, emax] in model from spectral
    // component (units: ph/cm2/s)
    double flux = spectral->flux(emin, emax);

    // Free spectral model if required
    if (free_spectral) delete spectral;

    // Return model flux
    return flux;
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

        // Get pointer on event list (circumvent const correctness)
        GCTAEventList* events =
            static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));

        // Loop over all models
        for (int i = 0; i < models.size(); ++i) {

            // Get data model (NULL if not a data model)
            const GModelData* model =
                dynamic_cast<const GModelData*>(models[i]);

            // If we have a data model that applies to the present observation
            // then simulate events
            if (model != NULL &&
                model->is_valid(obs->instrument(), obs->id())) {

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


/***********************************************************************//**
 * @brief Save event list in FITS format.
 *
 * Save the event list as a FITS file. The filename of the FITS file is
 * specified by the outfile parameter.
 ***************************************************************************/
void ctobssim::save_fits(void)
{
    // Save only if event list has not yet been saved and disposed
    if (!m_save_and_dispose) {

        // Get output filename
        m_outevents = (*this)["outevents"].filename();

        // Get CTA observation from observation container
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

        // Save observation into FITS file
        obs->save(m_outevents, clobber());

    } // endif: event list has not yet been saved and disposed

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
    m_outevents = (*this)["outevents"].filename();
    m_prefix  = (*this)["prefix"].string();

    // Issue warning if output filename has no .xml suffix
    std::string suffix = gammalib::tolower(m_outevents.substr(m_outevents.length()-4,4));
    if (suffix != ".xml") {
        log << "*** WARNING: Name of observation definition output file \""+
               m_outevents+"\"" << std::endl;
        log << "*** WARNING: does not terminate with \".xml\"." << std::endl;
        log << "*** WARNING: This is not an error, but might be misleading."
               " It is recommended" << std::endl;
        log << "*** WARNING: to use the suffix \".xml\" for observation"
               " definition files." << std::endl;
    }

    // Save only if event lists have not yet been saved and disposed
    if (!m_save_and_dispose) {

        // Loop over all observation in the container
        for (int i = 0; i < m_obs.size(); ++i) {

            // Get CTA observation
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

            // Handle only CTA observations
            if (obs != NULL) {

                // Continue only if there is an event list (it may have been disposed)
                if (obs->events()->size() != 0) {

                    // Set event output file name
                    std::string outfile = m_prefix + gammalib::str(i) + ".fits";

                    // Store output file name in observation
                    obs->eventfile(outfile);

                    // Save observation into FITS file
                    obs->save(outfile, clobber());

                }

            } // endif: observation was a CTA observations

        } // endfor: looped over observations

    } // endif: event list has not yet been saved and disposed

    // Save observations in XML file
    m_obs.save(m_outevents);

    // Return
    return;
}
