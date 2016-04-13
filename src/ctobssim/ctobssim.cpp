/***************************************************************************
 *                  ctobssim - Observation simulator tool                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2016 by Juergen Knoedlseder                         *
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
#include <typeinfo> 
#include "ctobssim.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                         "ctobssim::get_parameters()"
#define G_SIMULATE_SOURCE      "ctobssim::simulate_source(GCTAObservation*, "\
                                                    "GModels&, GRan&, GLog*)"
#define G_SIMULATE_INTERVAL    "ctobssim::simulate_interval(GCTAObservation*"\
                       ", GCTAEventList*, GCTAResponseIrf*, GModels&, GTime&"\
                            ", GTime&, GEnergy&, GEnergy&, GSkyDir&, double&"\
                                             ", double&, GRan&, GLog*, int&)"
#define G_GET_AREA      "ctobssim::get_area(GCTAObservation* obs, GEnergy&, "\
                                                                  "GEnergy&)"

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
 *
 * Constructs an empty ctobssim tool.
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
 * Constructs ctobssim tool from an observation container.
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
 * Constructs ctobssim tool using command line arguments for user parameter
 * setting.
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
 *
 * Constructs ctobssim tool from another ctobssim instance.
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
 *
 * Destructs ctobssim tool.
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
 * @param[in] app ctobssim tool.
 * @return ctobssim tool.
 *
 * Assigns ctobssim tool.
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
 * @brief Clear ctobssim tool
 *
 * Clears ctobssim tool.
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
 * @brief Run the ctobssim tool.
 *
 * Gets the user parameters, loops over all CTA observations in the
 * observation container, and simulate events for each observation. 
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

    // Determine the number of valid CTA observations, set energy dispersion
    // flag for all CTA observations and save old values in save_edisp vector
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
        log.header1(gammalib::number("Observation", m_obs.size()));
        log << m_obs << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1(gammalib::number("Simulate observation", m_obs.size()));
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

            // Write header for observation
            if (logTerse()) {
                wrklog.header3(get_obs_header(m_obs[i]));
            }

            // Get pointer on CTA observation
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

            // Skip observation if it's not CTA
            if (obs == NULL) {
                if (logTerse()) {
                    wrklog << " Skipping ";
                    wrklog << m_obs[i]->instrument();
                    wrklog << " observation" << std::endl;
                }
                continue;
            }

            // Skip observation if we have a binned observation
            if (obs->eventtype() == "CountsCube") {
                if (logTerse()) {
                    wrklog << " Skipping binned ";
                    wrklog << obs->instrument();
                    wrklog << " observation" << std::endl;
                }
                continue;
            }

            // Remove now all events from the event list but keep the
            // event list information such as ROI, Good Time Intervals,
            // energy boundaries. This will also keep additional columns
            // in an event list file.
            GCTAEventList* events = static_cast<GCTAEventList*>(obs->events());
            events->remove(0, events->size());

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
            if (logTerse()) {
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
                #pragma omp critical(ctobssim_run)
                {
                    //obs_clone.save(outfile, clobber());
                    obs->save(outfile, clobber());
                }

                // Dispose events
                obs->dispose_events();

            } // endif: save and dispose requested

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
    m_eslices     = 10;
    m_apply_edisp = false;
    m_max_rate    = 1.0e6;

    // Initialise protected members
    m_rans.clear();
    m_obs.clear();
    m_save_and_dispose = false;

    // Set fixed parameters
    m_max_photons = 1000000;  //!< Maximum number of photons / time slice

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
    m_eslices     = app.m_eslices;
    m_apply_edisp = app.m_apply_edisp;
    m_max_rate    = app.m_max_rate;

    // Copy protected members
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

        // Throw exception if counts cube is given
        require_inobs_nocube(G_GET_PARAMETERS);

        // Get observation container
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
    m_eslices     = (*this)["eslices"].integer();
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
 * @param[in,out] ran Random number generator.
 * @param[in] wrklog Pointer to logger.
 *
 * Simulate source events from a photon list for a given CTA observation and
 * all source models. The events are stored in form of an event list in the
 * observation.
 *
 * The method will loop over all Good Time Intervals (GTI) to perform the
 * simulations. Within a given GTI the requested energy range is split into
 * a number of energy slices so that the simulation area can be adapted to
 * the effective area of the instrument (this avoids simulating a large
 * number of sources photons at low energies while accepting only few of
 * them do to the small effective area). In case that the flux of a source
 * is large, the simulation may be done within time slices so that the
 * memory requirements won't get too large.
 *
 * This method does nothing if the observation pointer is NULL. It verifies
 * if the observation has a CTA IRF response.
 ***************************************************************************/
void ctobssim::simulate_source(GCTAObservation* obs,
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

        // Get pointer on event list
        GCTAEventList* events = static_cast<GCTAEventList*>(obs->events());

        // Get CTA response
        const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(obs->response());
        if (rsp == NULL) {
            std::string cls = std::string(typeid(obs->response()).name());
            std::string msg = "Response of type \""+cls+"\" is not a CTA "
                              "IRF response. Please make sure that a CTA "
                              "IRF response  is contained in the CTA "
                              "observation.";
            throw GException::invalid_value(G_SIMULATE_SOURCE, msg);
        }

        // Extract simulation region from event list ROI
        GSkyDir dir = events->roi().centre().dir();
        double  rad = events->roi().radius() + g_roi_margin;

        // Determine energy boundaries for simulation
        GEbounds ebounds = get_ebounds(events->ebounds());

        // Log simulation cone information
        if (logTerse()) {
            *wrklog << gammalib::parformat("Simulation cone");
            *wrklog << "RA=" << dir.ra_deg() << " deg";
            *wrklog << ", Dec=" << dir.dec_deg() << " deg";
            *wrklog << ", radius=" << rad << " deg" << std::endl;
        }

        // Initialise indentation for logging
        int indent = 0;

        // Initialise photon and event counters
        std::vector<int> nphotons(models.size(),0);
        std::vector<int> nevents(models.size(),0);

        // Loop over all Good Time Intervals
        for (int it = 0; it < events->gti().size(); ++it) {

            // Extract time interval
            GTime tmin = events->gti().tstart(it);
            GTime tmax = events->gti().tstop(it);

            // Log time interval. Increment indentation if there are
            // several Good Time Intervals.
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
            for (int ie = 0; ie <  ebounds.size(); ++ie) {

                // Set energy interval
                GEnergy emin = ebounds.emin(ie);
                GEnergy emax = ebounds.emax(ie);

                // Set true photon energy limits for simulation. If the
                // observation has energy dispersion then add a margin.
                GEnergy e_true_min = emin;
                GEnergy e_true_max = emax;
                if (rsp->use_edisp()) {
                    e_true_min = rsp->ebounds(e_true_min).emin();
                    e_true_max = rsp->ebounds(e_true_max).emax();
                }

                // Determine simulation area
                double area = get_area(obs, emin, emax);

                // Log energy range and simulation area
                if (logNormal()) {
                    *wrklog << gammalib::parformat("Photon energy range", indent);
                    *wrklog << e_true_min << " - " << e_true_max << std::endl;
                    *wrklog << gammalib::parformat("Event energy range", indent);
                    *wrklog << emin << " - " << emax << std::endl;
                }

                // Increment indentation if there are several energy
                // boundaries
                if (logNormal()) {
                    if (ebounds.size() > 1) {
                        indent++;
                        wrklog->indent(indent);
                    }
                }

                // Log simulation area
                if (logNormal()) {
                    *wrklog << gammalib::parformat("Simulation area", indent);
                    *wrklog << area << " cm2" << std::endl;
                }

                // Save state of event counter before doing the simulation
                int nevents_before = events->size();

                // Simulate events for this time and energy interval
                simulate_interval(obs, rsp, events, models, tmin, tmax,
                                  e_true_min, e_true_max, dir, rad, area,
                                  ran, wrklog, indent, nphotons, nevents);

                // Log simulation results
                if (logNormal()) {
                    *wrklog << gammalib::parformat("MC source events", indent);
                    *wrklog << events->size() - nevents_before;
                    *wrklog << " (all source models)";
                    *wrklog << std::endl;
                }

                // Reset indentation if there were several energy boundaries
                if (logNormal()) {
                    if (ebounds.size() > 1) {
                        indent--;
                        wrklog->indent(indent);

                    }
                }

            } // endfor: looped over all energy boundaries

            // Reset indentation if there were several Good Time Intervals
            if (logNormal()) {
                if (events->gti().size() > 1) {
                    indent--;
                    wrklog->indent(indent);
                }
            }

        } // endfor: looped over all time intervals

        // Reset indentation
        wrklog->indent(0);

        // Log simulation summary
        if (logTerse()) {
            for (int i = 0; i < models.size(); ++i) {
                const GModelSky* model = dynamic_cast<const GModelSky*>(models[i]);
                if (model == NULL) {
                    continue;
                }
                if (!model->is_valid(obs->instrument(), obs->id())) {
                    continue;
                }
                *wrklog << gammalib::parformat("MC source photons");
                *wrklog << nphotons[i];
                if (model->name().length() > 0) {
                    *wrklog << " [" << model->name() << "]";
                }
                *wrklog << std::endl;
                *wrklog << gammalib::parformat("MC source events");
                *wrklog << nevents[i];
                if (model->name().length() > 0) {
                    *wrklog << " [" << model->name() << "]";
                }
                *wrklog << std::endl;
            }
        }

    } // endif: observation pointer was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Simulate source events for a time and energy interval
 *
 * @param[in] obs Pointer on CTA observation.
 * @param[in] rsp Pointer on CTA IRF response.
 * @param[in,out] events Pointer on CTA event list.
 * @param[in] models Model list.
 * @param[in] tmin Start time.
 * @param[in] tmax Stop time.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 * @param[in] dir Simulation cone centre.
 * @param[in] rad Simulation cone radius (degrees).
 * @param[in] area Simulation area (cm^2).
 * @param[in,out] ran Random number generator.
 * @param[in] wrklog Pointer to logger.
 * @param[in,out] indent Logger indent.
 * @param[in,out] nphotons Number of photons for all models.
 * @param[in,out] nevents Number of events for all models.
 *
 * Simulate source events for a time and energy interval.
 ***************************************************************************/
void ctobssim::simulate_interval(GCTAObservation*       obs,
                                 const GCTAResponseIrf* rsp,
                                 GCTAEventList*         events,
                                 const GModels&         models,
                                 const GTime&           tmin,
                                 const GTime&           tmax,
                                 const GEnergy&         emin,
                                 const GEnergy&         emax,
                                 const GSkyDir&         dir,
                                 const double&          rad,
                                 const double&          area,
                                 GRan&                  ran,
                                 GLog*                  wrklog,
                                 int&                   indent,
                                 std::vector<int>&      nphotons,
                                 std::vector<int>&      nevents)
{
    // Loop over all models
    for (int i = 0; i < models.size(); ++i) {

        // Get sky model (NULL if not a sky model)
        const GModelSky* model = dynamic_cast<const GModelSky*>(models[i]);

        // If the model is not a sky model then skip the model
        if (model == NULL) {
            continue;
        }

        // If the model does not apply to the instrument and observation
        // identifier then skip the model
        if (!model->is_valid(obs->instrument(), obs->id())) {
            continue;
        }

        // Determine duration of a time slice by limiting the number of
        // simulated photons to m_max_photons. The photon rate is estimated
        // from the model flux and used to set the duration of the time
        // slice.
        double flux     = get_model_flux(model, emin, emax, dir, rad,
                                         indent, wrklog);
        double rate     = flux * area;
        double duration = 1800.0;           // default: 1800 sec
        if (rate > 0.0) {
            duration = m_max_photons / rate;
            if (duration < 1.0) {           // not <1 sec
                duration = 1.0;
            }
            else if (duration > 180000.0) { // not >50 hr
                duration = 180000.0;
            }
        }

        // Skip model if photon rate is 0
        if (rate <= 0.0) {
            continue;
        }

        // Log photon rate
        if (logNormal()) {
            *wrklog << gammalib::parformat("Photon rate", indent);
            *wrklog << rate << " photons/s";
            if (model->name().length() > 0) {
                *wrklog << " [" << model->name() << "]";
            }
            *wrklog << std::endl;
        }

        // If photon rate exceeds the maximum photon rate that is allowed
        // then throw an exception
        if (rate > m_max_rate) {
            std::string mod = (model->name().length() > 0) ?
                               model->name() : "Unknown";
            std::string msg = "Photon rate "+gammalib::str(rate)+
                              " photons/s for model \""+mod+"\" exceeds "
                              "maximum allowed photon rate of "+
                              gammalib::str(m_max_rate)+" photons/s. "
                              "Please check the parameters of model "
                              "\""+mod+"\" or increase the value of the "
                              "\"maxrate\" parameter.";
            throw GException::invalid_value(G_SIMULATE_INTERVAL, msg);
        }

        // To reduce memory requirements we split long time intervals into
        // several time slices
        GTime tstart = tmin;
        GTime tstop  = tstart + duration;

        // Save state of photon and event counters before doing the
        // simulation
        int nphotons_before = nphotons[i];
        int nevents_before  = nevents[i];

        // Loop over time slices
        while (tstart < tmax) {

            // Make sure that tstop <= tmax
            if (tstop > tmax) {
                tstop = tmax;
            }

            // Log time slice
            if (logExplicit()) {
                if (tmax - tmin > duration) {
                    indent++;
                    wrklog->indent(indent);
                }
                *wrklog << gammalib::parformat("Time slice", indent);
                *wrklog << tstart.convert(m_cta_ref) << " - ";
                *wrklog << tstop.convert(m_cta_ref) << " s";
                if (model->name().length() > 0) {
                    *wrklog << " [" << model->name() << "]";
                }
                *wrklog << std::endl;
            }

            // Simulate time slice
            simulate_time_slice(obs, rsp, events, model, tstart, tstop,
                                emin, emax, dir, rad, area,
                                ran, wrklog, indent,
                                nphotons[i], nevents[i]);

            // Go to next time slice
            tstart = tstop;
            tstop  = tstart + duration;

            // Reset indentation
            if (logExplicit()) {
                if (tmax - tmin > duration) {
                    indent--;
                    wrklog->indent(indent);
                }
            }

        } // endwhile: looped over time slices

        // Log simulation results
        if (logNormal()) {
            *wrklog << gammalib::parformat("MC source photons", indent);
            *wrklog << nphotons[i] - nphotons_before;
            if (model->name().length() > 0) {
                *wrklog << " [" << model->name() << "]";
            }
            *wrklog << std::endl;
            *wrklog << gammalib::parformat("MC source events", indent);
            *wrklog << nevents[i] - nevents_before;
            if (model->name().length() > 0) {
                *wrklog << " [" << model->name() << "]";
            }
            *wrklog << std::endl;

        }
            
    } // endfor: looped over models

    // Return
    return;
}


/***********************************************************************//**
 * @brief Simulate source events for a time slice
 *
 * @param[in] obs Pointer on CTA observation.
 * @param[in,out] events Pointer on CTA event list.
 * @param[in] models Model list.
 * @param[in] tmin Start time.
 * @param[in] tmax Stop time.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 * @param[in] dir Simulation cone centre.
 * @param[in] rad Simulation cone radius (degrees).
 * @param[in] area Simulation area (cm^2).
 * @param[in,out] ran Random number generator.
 * @param[in] wrklog Pointer to logger.
 * @param[in,out] indent Logger indent.
 * @param[in,out] nphotons Number of photons.
 * @param[in,out] nevents Number of events.
 *
 * Simulate source events for a time slice.
 ***************************************************************************/
void ctobssim::simulate_time_slice(GCTAObservation*       obs,
                                   const GCTAResponseIrf* rsp,
                                   GCTAEventList*         events,
                                   const GModelSky*       model,
                                   const GTime&           tstart,
                                   const GTime&           tstop,
                                   const GEnergy&         emin,
                                   const GEnergy&         emax,
                                   const GSkyDir&         dir,
                                   const double&          rad,
                                   const double&          area,
                                   GRan&                  ran,
                                   GLog*                  wrklog,
                                   int&                   indent,
                                   int&                   nphotons,
                                   int&                   nevents)
{
    // Get photons
    GPhotons photons = model->mc(area, dir, rad, emin, emax, tstart, tstop,
                                 ran);

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

        // Simulate event. Note that this method includes the deadtime
        // correction.
        GCTAEventAtom* event = rsp->mc(area, photons[i], *obs, ran);

        // Use event only if it exists and if it falls within ROI, the
        // energy interval and the time slice
        if (event != NULL) {
            if (events->roi().contains(*event) &&
                event->energy() >= emin &&
                event->energy() <= emax &&
                event->time() >= tstart &&
                event->time() <= tstop) {
                event->event_id(m_event_id);
                events->append(*event);
                m_event_id++;
                nevents++;
            }
            delete event;
        }

    } // endfor: looped over events

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get energy boundaries
 *
 * @param[in] ebounds Energy boundaries of events.
 * @return Energy boundaries for simulation.
 *
 * Get the energy boundaries for the source simulation.
 *
 * @todo Need to take care of the (rare) case that there are distinct
 *       energy boundaries in the event list. In that case we need to
 *       create mixed energy boundaries.
 ***************************************************************************/
GEbounds ctobssim::get_ebounds(const GEbounds& ebounds) const
{
    // Set energy bins
    GEbounds ebins(m_eslices, ebounds.emin(), ebounds.emax());
 
    // Return energy bins
    return (ebins);
}


/***********************************************************************//**
 * @brief Get simulation area (cm^2)
 *
 * @param[in] obs Pointer on CTA observation.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 * @return Simulation area (cm^2).
 *
 * Get the simulation area for an energy interval in units of cm^2.
 ***************************************************************************/
double ctobssim::get_area(GCTAObservation* obs,
                          const GEnergy&   emin,
                          const GEnergy&   emax) const
{
    // Get CTA response
    const GCTAResponseIrf* rsp =
                      dynamic_cast<const GCTAResponseIrf*>(obs->response());
    if (rsp == NULL) {
        std::string cls = std::string(typeid(rsp).name());
        std::string msg = "Response of type \""+cls+"\" is not a CTA IRF "
                          "response. Please make sure that a CTA IRF "
                          "response  is contained in the CTA observation.";
        throw GException::invalid_value(G_GET_AREA, msg);
    }

    // Compute effective area at minimum and maximum energy
    double area_emin = rsp->aeff()->max(emin.log10TeV(), 0.0, 0.0);
    double area_emax = rsp->aeff()->max(emax.log10TeV(), 0.0, 0.0);

    // Use maximum of both
    double area = (area_emin > area_emax) ? area_emin : area_emax;

    // Multiply by security factor
    area *= 2.0;

    // Return simulation area
    return area;
}


/***********************************************************************//**
 * @brief Determine sky model flux
 *
 * @param[in] model Sky model.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 * @param[in] centre Centre of region for photon rate determination.
 * @param[in] radius Radius of region for photon rate determination (degrees).
 * @param[in] indent Indent for logging.
 * @param[in,out] wrklog Pointer to logger.
 * @return Model flux (photons/cm2/sec).
 ***************************************************************************/
double ctobssim::get_model_flux(const GModelSky* model,
                                const GEnergy&   emin,
                                const GEnergy&   emax,
                                const GSkyDir&   centre,
                                const double&    radius,
                                const int&       indent,
                                GLog*            wrklog)
{
    // Initialise flux
    double flux = 0.0;

    // Determine the spatial model normalization within the simulation
    // cone and check whether the model will produce any photons in that
    // cone.
    double norm      = model->spatial()->mc_norm(centre, radius);
    bool   use_model = (norm > 0.0) ? true : false;
    if (logNormal()) {
        if (use_model) {
            *wrklog << gammalib::parformat("Use model", indent);
        }
        else {
            *wrklog << gammalib::parformat("Skip model", indent);
        }
        if (model->name().length() > 0) {
            *wrklog << model->name();
        }
        *wrklog << std::endl;
        *wrklog << gammalib::parformat("Normalization", indent);
        *wrklog << norm;
        if (model->name().length() > 0) {
            *wrklog << " [" << model->name() << "]";
        }
        *wrklog << std::endl;
    }

    // Continue only if model overlaps with simulation region
    if (use_model) {

        // Initialise de-allocation flag
        bool free_spectral = false;

        // Get pointer to spectral model
        const GModelSpectral* spectral = model->spectral();

        // If the spatial model is a diffuse cube then create a node
        // function spectral model that is the product of the diffuse
        // cube node function and the spectral model evaluated at the
        // energies of the node function
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
                double  value     = spectral->eval(energy, time);
                nodes->intensity(i, value * intensity);
            }

            // Signal that node function needs to be de-allocated later
            free_spectral = true;

            // Set the spectral model pointer to the node function
            spectral = nodes;

        } // endif: spatial model was a diffuse cube
    
        // Compute flux within [emin, emax] in model from spectral
        // component (units: ph/cm2/s)
        double flux0 = spectral->flux(emin, emax);
        flux         = flux0 * norm;

        // Dump flux
        if (logNormal()) {
            *wrklog << gammalib::parformat("Flux", indent);
            *wrklog << flux0;
            if (model->name().length() > 0) {
                *wrklog << " [" << model->name() << "]";
            }
            *wrklog << " photons/cm2/s" << std::endl;
            *wrklog << gammalib::parformat("Normalized flux", indent);
            *wrklog << flux;
            if (model->name().length() > 0) {
                *wrklog << " [" << model->name() << "]";
            }
            *wrklog << " photons/cm2/s" << std::endl;
        }

        // Free spectral model if required
        if (free_spectral) delete spectral;

    } // endif: model overlaps with simulation region

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

                    // Initialise statistics
                    int n_appended    = 0;
                    int n_outside_roi = 0;
                    
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

                            // Increment number of appended events
                            n_appended++;

                        } // endif: event was within ROI

                        // ... otherwise increment outside ROI counter
                        else {
                            n_outside_roi++;
                        }

                    } // endfor: looped over all events

                    // Dump simulation results
                    if (logNormal()) {
                        *wrklog << gammalib::parformat("MC events outside ROI");
                        *wrklog << n_outside_roi << std::endl;
                        *wrklog << gammalib::parformat("MC background events");
                        *wrklog << n_appended << std::endl;
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

        // Log filename
        if (logTerse()) {
            log << gammalib::parformat("Event list file");
            log << m_outevents << std::endl;
        }

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
    m_prefix    = (*this)["prefix"].string();

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
                    std::string outfile = m_prefix + gammalib::str(i) +
                                          ".fits";

                    // Store output file name in observation
                    obs->eventfile(outfile);

                    // Log filename
                    if (logTerse()) {
                        log << gammalib::parformat("Event list file");
                        log << outfile << std::endl;
                    }

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
