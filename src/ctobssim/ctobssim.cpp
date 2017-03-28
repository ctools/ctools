/***************************************************************************
 *                  ctobssim - Observation simulator tool                  *
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
#include "GFits.hpp"


/* __ Method name definitions ____________________________________________ */
#define G_GET_PARAMETERS                         "ctobssim::get_parameters()"
#define G_SIMULATE_SOURCE      "ctobssim::simulate_source(GCTAObservation*, "\
                                                    "GModels&, GRan&, GLog*)"
#define G_SIMULATE_INTERVAL    "ctobssim::simulate_interval(GCTAObservation*"\
                       ", GCTAEventList*, GCTAResponseIrf*, GModels&, GTime&"\
                           ", GTime&, GEnergy&, GEnergy&, GEnergy&, GEnergy&"\
                          ", GSkyDir&, double&, double&, GRan&, GLog*, int&)"
#define G_GET_AREA      "ctobssim::get_area(GCTAObservation* obs, GEnergy&, "\
                                                                  "GEnergy&)"

/* __ Constants __________________________________________________________ */
const double g_roi_margin = 0.5;      //!< Simulation radius margin (degrees)

/* __ Debug definitions __________________________________________________ */
//#define G_SOURCE_DEBUG
//#define G_BACKGROUND_DEBUG

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
ctobssim::ctobssim(void) : ctobservation(CTOBSSIM_NAME, CTOBSSIM_VERSION)
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
          ctobservation(CTOBSSIM_NAME, CTOBSSIM_VERSION, obs)
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
 * Constructs ctobssim tool using command line arguments for user parameter
 * setting.
 ***************************************************************************/
ctobssim::ctobssim(int argc, char *argv[]) :
          ctobservation(CTOBSSIM_NAME, CTOBSSIM_VERSION, argc, argv)
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
ctobssim::ctobssim(const ctobssim& app) : ctobservation(app)
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
 * @brief Clear ctobssim tool
 *
 * Clears ctobssim tool.
 ***************************************************************************/
void ctobssim::clear(void)
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

    // Special mode: if read ahead is specified we know that we called
    // the execute() method, hence files are saved immediately and event
    // lists are disposed afterwards.
    if (read_ahead()) {
        m_save_and_dispose = true;
    }

    // Set energy dispersion flags of all CTA observations and save old
    // values in save_edisp vector
    std::vector<bool> save_edisp = set_edisp(m_obs, m_apply_edisp);

    // Determine the number of valid CTA observations
    int n_observations = 0;
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            n_observations++;
        }
    }

    // If more than a single observation has been handled then make sure that
    // an XML file will be used for storage
    if (n_observations > 1) {
        m_use_xml = true;
    }

    // Write execution model into logger
    log_header1(NORMAL, "Execution mode");
    std::string mode = (m_save_and_dispose)
                       ? "Save and dispose (reduces memory needs)"
                       : "Keep events in memory";
    std::string xml  = (m_use_xml)
                       ? "Write Observation Definition XML file"
                       : "Write single event list FITS file";
    log_value(NORMAL, "Event list management", mode);
    log_value(NORMAL, "Output format", xml);

    // Write seed values into logger
    log_header1(NORMAL, "Seed values");
    for (int i = 0; i < m_rans.size(); ++i) {
        log_value(NORMAL, "Seed "+gammalib::str(i),
                  gammalib::str(m_rans[i].seed()));
    }

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write header
    log_header1(TERSE, gammalib::number("Simulate observation", m_obs.size()));

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
        wrklog.buffer_size(10000000);

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

                // Set event output file name
                std::string outfile = this->outfile(i);

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

    // Restore energy dispersion flags of all CTA observations
    restore_edisp(m_obs, save_edisp);

    // Optionally publish event list(s)
    if ((*this)["publish"].boolean()) {
        publish();
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
    // Write header into logger
    log_header1(TERSE, gammalib::number("Save observation", m_obs.size()));

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
 * @brief Publish event lists
 *
 * @param[in] name Event list name.
 ***************************************************************************/
void ctobssim::publish(const std::string& name)
{
    // Write header into logger
    log_header1(TERSE, gammalib::number("Publish event list", m_obs.size()));

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Handle only CTA observations
        if (obs != NULL) {

            // Continue only if there is an event list
            if (obs->events()->size() != 0) {

                // Set default name if user name is empty
                std::string user_name(name);
                if (user_name.empty()) {
                    user_name = CTOBSSIM_NAME;
                }

                // If there are several event lists then add an index
                if (m_use_xml) {
                    user_name += gammalib::str(i);
                }

                // Write event list name into logger
                log_value(NORMAL, "Event list name", user_name);

                // Write events into in-memory FITS file
                GFits fits;
                obs->write(fits);

                // Publish
                fits.publish(gammalib::extname_cta_events, user_name);

            } // endif: there were events

        } // endif: observation was a CTA observation

    } // endfor: looped over observations

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
    m_startindex  = 1;
    m_seed        = 1;
    m_eslices     = 10;
    m_apply_edisp = false;
    m_max_rate    = 1.0e6;

    // Initialise protected members
    m_rans.clear();
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
    m_startindex  = app.m_startindex;
    m_seed        = app.m_seed;
    m_eslices     = app.m_eslices;
    m_apply_edisp = app.m_apply_edisp;
    m_max_rate    = app.m_max_rate;

    // Copy protected members
    m_max_photons      = app.m_max_photons;
    m_rans             = app.m_rans;
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

    // ... otherwise add response information and energy boundaries in case
    // that they are missing
    else {
        setup_observations(m_obs);
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
        m_outevents  = (*this)["outevents"].filename();
        m_prefix     = (*this)["prefix"].string();
        m_startindex = (*this)["startindex"].integer();
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

    // Write parameters into logger
    log_parameters(TERSE);

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
    // Debug code: signal that we step into the simulate_source method
    #if defined(G_SOURCE_DEBUG)
    std::cout << "ctobssim::simulate_source: in" << std::endl;
    #endif

    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // If no logger is specified then use the default logger
        if (wrklog == NULL) {
            wrklog = &log;
        }

        // Get pointer on event list
        GCTAEventList* events = static_cast<GCTAEventList*>(obs->events());

        // Get CTA response
        const GCTAResponseIrf* rsp =
              dynamic_cast<const GCTAResponseIrf*>(obs->response());
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

                // Set reconstructed energy interval
                GEnergy ereco_min = ebounds.emin(ie);
                GEnergy ereco_max = ebounds.emax(ie);

                // Set true photon energy limits for simulation. If the
                // observation has energy dispersion then add a margin.
                GEnergy etrue_min = ereco_min;
                GEnergy etrue_max = ereco_max;
                if (rsp->use_edisp()) {
                    etrue_min = rsp->ebounds(etrue_min).emin();
                    etrue_max = rsp->ebounds(etrue_max).emax();
                }

                // Determine simulation area
                double area = get_area(obs, etrue_min, etrue_max);

                // Log energy range and simulation area
                if (logNormal()) {
                    *wrklog << gammalib::parformat("Photon energy range", indent);
                    *wrklog << etrue_min << " - " << etrue_max << std::endl;
                    *wrklog << gammalib::parformat("Event energy range", indent);
                    *wrklog << ereco_min << " - " << ereco_max << std::endl;
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
                                  etrue_min, etrue_max, ereco_min, ereco_max,
                                  dir, rad, area,
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

    // Debug code: signal that we are down with the simulate_source method
    #if defined(G_SOURCE_DEBUG)
    std::cout << "ctobssim::simulate_source: out" << std::endl;
    #endif

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
 * @param[in] etrue_min Minimum true energy.
 * @param[in] etrue_max Maximum true energy.
 * @param[in] ereco_min Minimum reconstructed energy.
 * @param[in] ereco_max Maximum reconstructed energy.
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
                                 const GEnergy&         etrue_min,
                                 const GEnergy&         etrue_max,
                                 const GEnergy&         ereco_min,
                                 const GEnergy&         ereco_max,
                                 const GSkyDir&         dir,
                                 const double&          rad,
                                 const double&          area,
                                 GRan&                  ran,
                                 GLog*                  wrklog,
                                 int&                   indent,
                                 std::vector<int>&      nphotons,
                                 std::vector<int>&      nevents)
{
    // Debug code: signal that we step into the simulate_interval method
    #if defined(G_SOURCE_DEBUG)
    std::cout << "ctobssim::simulate_interval: in";
    std::cout << " Etrue=[" << etrue_min << "," << etrue_max << "]";
    std::cout << " Ereco=[" << ereco_min << "," << ereco_max << "]";
    std::cout << std::endl;
    #endif

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
        double flux     = get_model_flux(model, etrue_min, etrue_max, dir, rad,
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
                                etrue_min, etrue_max, ereco_min, ereco_max,
                                dir, rad, area,
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

    // Debug code: signal that we are down with the simulate_interval method
    #if defined(G_SOURCE_DEBUG)
    std::cout << "ctobssim::simulate_interval: out" << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Simulate source events for a time slice
 *
 * @param[in] obs Pointer on CTA observation.
 * @param[in] rsp Pointer on CTA response.
 * @param[in,out] events Pointer on CTA event list.
 * @param[in] model Model.
 * @param[in] tstart Start time.
 * @param[in] tstop Stop time.
 * @param[in] etrue_min Minimum true energy.
 * @param[in] etrue_max Maximum true energy.
 * @param[in] ereco_min Minimum reconstructed energy.
 * @param[in] ereco_max Maximum reconstructed energy.
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
                                   const GEnergy&         etrue_min,
                                   const GEnergy&         etrue_max,
                                   const GEnergy&         ereco_min,
                                   const GEnergy&         ereco_max,
                                   const GSkyDir&         dir,
                                   const double&          rad,
                                   const double&          area,
                                   GRan&                  ran,
                                   GLog*                  wrklog,
                                   int&                   indent,
                                   int&                   nphotons,
                                   int&                   nevents)
{
    // Debug code: signal that we step into the model MC method
    #if defined(G_SOURCE_DEBUG)
    std::cout << "ctobssim::simulate_time_slice: model->mc in" << std::endl;
    #endif

    // Get photons
    GPhotons photons = model->mc(area, dir, rad, etrue_min, etrue_max,
                                 tstart, tstop, ran);

    // Debug code: signal that we came back
    #if defined(G_SOURCE_DEBUG)
    std::cout << "ctobssim::simulate_time_slice: model->mc out" << std::endl;
    #endif

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

        // Debug code: signal that we step into the response MC method
        #if defined(G_SOURCE_DEBUG)
        std::cout << "ctobssim::simulate_time_slice: rsp->mc in" << std::endl;
        #endif

        // Simulate event. Note that this method includes the deadtime
        // correction.
        GCTAEventAtom* event = rsp->mc(area, photons[i], *obs, ran);

        // Debug code: signal that we came back
        #if defined(G_SOURCE_DEBUG)
        std::cout << "ctobssim::simulate_time_slice: rsp->mc out" << std::endl;
        #endif

        // Use event only if it exists and if it falls within ROI, the
        // reconstructed energy interval and the time slice
        if (event != NULL) {
            if (events->roi().contains(*event) &&
                event->energy() >= ereco_min &&
                event->energy() <= ereco_max &&
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
 * @param[in] emin Minimum true energy.
 * @param[in] emax Maximum true energy.
 * @return Simulation area (cm^2).
 *
 * Get the simulation area for an energy interval in units of cm^2. This is
 * done by extracting the maximum effective area value within the energy
 * range [emin,emax] and by multiplying this value by 2 for security. The
 * effective area is sampled at 100 energy values within the energy interval.
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
    const int nbins   = 10;
    double    logE    = emin.log10TeV();
    double    logEbin = (emax.log10TeV() - logE)/double(nbins-1);
    double    area    = 0.0;
    for (int i = 0; i < nbins; ++i, logE += logEbin) {
        double aeff = rsp->aeff()->max(logE, 0.0, 0.0);
        if (aeff > area) {
            area = aeff;
        }
    }

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
    // Debug code: signal that we step into the get_model_flux method
    #if defined(G_SOURCE_DEBUG)
    std::cout << "ctobssim::get_model_flux: in";
    std::cout << " Etrue=[" << emin << "," << emax << "]";
    std::cout << std::endl;
    #endif

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
        GModelSpatialDiffuseCube* cube =
                     dynamic_cast<GModelSpatialDiffuseCube*>(model->spatial());
        if (cube != NULL) {

            // Set MC cone
            cube->set_mc_cone(centre, radius);

            // Allocate node function to replace the spectral component
            GModelSpectralNodes* nodes =
                          new GModelSpectralNodes(cube->spectrum());
            for (int i = 0; i < nodes->nodes(); ++i) {
                GEnergy energy    = nodes->energy(i);
                double  intensity = nodes->intensity(i);
                double  value     = spectral->eval(energy);
                nodes->intensity(i, value * intensity);
            }

            // Signal that node function needs to be de-allocated later
            free_spectral = true;

            // Set the spectral model pointer to the node function
            spectral = nodes;

            // Kluge: if there are no nodes then the spectral->flux method
            // will throw an exception. We therefore set here the use_model
            // flag to false in case that there are no spectral nodes
            if (nodes->nodes() == 0) {
                use_model = false;
            }

        } // endif: spatial model was a diffuse cube

        // Compute flux within [emin, emax] in model from spectral
        // component (units: ph/cm2/s)
        double flux0 = (use_model) ? spectral->flux(emin, emax) : 0.0;
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

    // Debug code: signal that we step out of the get_model_flux method
    #if defined(G_SOURCE_DEBUG)
    std::cout << "ctobssim::get_model_flux: in" << std::endl;
    #endif

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

                // Debug code: signal that we step into the response MC method
                #if defined(G_BACKGROUND_DEBUG)
                std::cout << "ctobssim::simulate_background: model->mc in" << std::endl;
                #endif

                // Get simulated CTA event list. Note that this method
                // includes the deadtime correction.
                GCTAEventList* list =
                     dynamic_cast<GCTAEventList*>(model->mc(*obs, ran));

                // Debug code: signal that we came back
                #if defined(G_BACKGROUND_DEBUG)
                std::cout << "ctobssim::simulate_background: model->mc out" << std::endl;
                #endif

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
    // Save only if event list has not yet been saved and disposed and if
    // there are observations
    if (!m_save_and_dispose && m_obs.size() > 0) {

        // Get output filename
        m_outevents = (*this)["outevents"].filename();

        // Get CTA observation from observation container
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

        // Log filename
        log_value(NORMAL, "Event list file", m_outevents);

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
    // Get output filename, prefix and start index
    m_outevents  = (*this)["outevents"].filename();

    // Issue warning if output filename has no .xml suffix
    log_string(TERSE, warn_xml_suffix(m_outevents));

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
                    std::string outfile = this->outfile(i);

                    // Store output file name in observation
                    obs->eventfile(outfile);

                    // Log filename
                    log_value(NORMAL, "Event list file", outfile);

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


/***********************************************************************//**
 * @brief Return output filename
 *
 * @param[in] index Observation index
 * @return Output filename
 *
 * Return output filename for observation with @p index.
 ***************************************************************************/
std::string ctobssim::outfile(const int& index)
{
    // Initialise output filename
    std::string outfile;

    // If multiple observations are handled then build the filename from
    // prefix and observation index plus startindex. The format is
    // [prefix]NNNNNN.fits, where NNNNNN is a 6-digit integer
    if (m_use_xml) {

        // Get prefix and start index
        m_prefix     = (*this)["prefix"].string();
        m_startindex = (*this)["startindex"].integer();

        // Build filename
        char buffer[256];
        std::sprintf(buffer, "%s%6.6d.fits", m_prefix.c_str(),
                                             index + m_startindex);

        // Set output filename
        outfile = std::string(buffer);

    }

    // ... otherwise use the outfile parameter
    else {

        // Get output event list file name
        m_outevents = (*this)["outevents"].filename();

        // Set output filename
        outfile = std::string(m_outevents);
    }

    // Return
    return outfile;
}
