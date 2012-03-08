/***************************************************************************
 *                ctobssim - CTA observation simulator tool                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2012 by Juergen Knoedlseder                         *
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
 * @author J. Knoedlseder
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
 * This constructor creates an instance of the class that is initialised from
 * an observation container.
 ***************************************************************************/
ctobssim::ctobssim(GObservations obs) : GApplication(CTOBSSIM_NAME, CTOBSSIM_VERSION)
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
 * @brief Clear instance
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
    // Read ahead output filename so that it gets dumped correctly in the
    // parameters log
    m_outfile = (*this)["outfile"].filename();

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
    if (logDebug())
        log.cout(true);

    // Get parameters
    get_parameters();

    // Write input parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
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

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(&m_obs[i]);

        // Continue only if observation is a CTA observation
        if (obs != NULL) {

            // Write header for observation
            if (logTerse()) {
                if (obs->name().length() > 1) {
                    log.header3("Observation "+obs->name());
                }
                else {
                    log.header3("Observation");
                }
            }

            // Simulate source events
            simulate_source(obs, m_obs.models());

            // Simulate source events
            simulate_background(obs, m_obs.models());

        } // endif: CTA observation found

    } // endfor: looped over observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save simulated observation
 *
 * This method saves the results.
 *
 * @todo Implement XML file output in case that on observation list is
 *       provided on input.
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

    // Get output filename
    m_outfile = (*this)["outfile"].filename();

    // Loop over all observation in the container
    int file_num = 0;
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(&m_obs[i]);

        // Save only if observation is a CTA observation
        if (obs != NULL) {

            // Set filename. If more than one file will be created an
            // index "_xxx" will be appended.
            std::string filename = m_outfile;
            if (file_num > 0) {
                filename += "_"+str(file_num);
            }

            // Save file
            obs->save(filename, clobber());

            // Increment file number
            file_num++;

        } // endif: observation was a CTA observation

    } // endfor: looped over files

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

    } // endif: there was no observation in the container

    // Get other parameters
    m_seed = (*this)["seed"].integer();

    // Initialise random number generator
    m_ran.seed(m_seed);

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
        GCTARoi     roi;
        GCTAInstDir instdir;
        instdir.radec_deg(m_ra, m_dec);
        roi.centre(instdir);
        roi.radius(m_rad);

        // Set GTI
        GGti  gti;
        GTime tstart;
        GTime tstop;
        tstart.met(m_tmin);
        tstop.met(m_tmax);
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
        obs->events(&events);
        
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
 * @param[in] photons Photon list.
 *
 * @exception GCTAException::no_pointing
 *            No valid pointing found in CTA observation
 * @exception GCTAException::no_response
 *            No valid response found in CTA observation
 *
 * Simulate source events from a photon list for a given CTA observation.
 * The events are stored in as event list in the observation.
 *
 * This method does nothing if the observation pointer is NULL. It also
 * verifies if the observation has a valid pointing and response.
 *
 * @todo Add margins to the data selections so that we generate also photons
 *       that fall slightly outside the ROI and the energy interval so that
 *       they can be scatter within the ROI and energy interval by the
 *       PSF and energy dispersion.
 ***************************************************************************/
void ctobssim::simulate_source(GCTAObservation* obs, const GModels& models)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Get pointer on CTA pointing. As we have no time dependent
        // pointings implement we just determine the pointing from the
        // beginning of the observation. Throw an exception if the pointing
        // is not defined.
        GCTAPointing* pnt = obs->pointing();
        if (pnt == NULL) {
            throw GCTAException::no_pointing(G_SIMULATE_SOURCE);
        }

        // Get pointer on CTA response. Throw an exception if the response
        // is not defined.
        GCTAResponse* rsp = obs->response();
        if (rsp == NULL) {
            throw GCTAException::no_response(G_SIMULATE_SOURCE);
        }

        // Make sure that the observation holds a CTA event list. If this
        // is not the case then allocate and attach a CTA event list now.
        if (dynamic_cast<const GCTAEventList*>(obs->events()) == NULL) {
            set_list(obs);
        }

        // Get pointer on event list (circumvent const correctness)
        GCTAEventList* events = (GCTAEventList*)(obs->events());

        // Extract ROI
        GSkyDir dir = events->roi().centre().skydir();
        double  rad = events->roi().radius();

        // Dump simulation cone information
        if (logNormal()) {
            log << parformat("Simulation area");
            log << str(m_area) << " cm2" << std::endl;
            log << parformat("Simulation cone");
            log << "RA=" << dir.ra_deg() << " deg";
            log << ", Dec=" << dir.dec_deg() << "deg";
            log << ", r=" << rad << " deg" << std::endl;
        }

        // Initialise indentation for logging
        int indent = 0;

        // Loop over all Good Time Intervals
        for (int it = 0; it <  events->gti().size(); ++it) {

            // Extract time interval
            GTime tmin = events->gti().tstart(it);
            GTime tmax = events->gti().tstop(it);
            
            // Dump time interval
            if (logNormal()) {
                if (events->gti().size() > 1) {
                    indent++;
                    log.indent(indent);
                }
                log << parformat("Time interval");
                log << tmin << " - " << tmax << std::endl;
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
                        log.indent(indent);
                    }
                    log << parformat("Energy range");
                    log << emin << " - " << emax << std::endl;
                }

                // Loop over all sky models
                for (int i = 0; i < models.size(); ++i) {
                
                    // Get sky model (NULL if not a sky model)
                    const GModelSky* model = dynamic_cast<const GModelSky*>(&models[i]);

                    // If we have a sky model then simulate photons
                    if (model != NULL) {

                        // To reduce memory requirements we split long time
                        // intervals into several slices. The maximum length
                        // of a slice is determined by the fixed parameter
                        // m_time_max
                        GTime tstart = tmin;
                        GTime tstop  = tstart + m_time_max;
                        
                        // Initialise cumulative photon counters
                        int nphotons = 0;
                        int ndeadc   = 0;
                        
                        // Loop over time slices
                        while (tstart < tmax) {
            
                            // Make sure that tstop <= tmax
                            if (tstop > tmax) {
                                tstop = tmax;
                            }

                            // Dump time slice
                            if (logExplicit()) {
                                if (tmax - tmin > m_time_max) {
                                    indent++;
                                    log.indent(indent);
                                }
                                log << parformat("Time slice");
                                log << tstart << " - " << tstop << std::endl;
                            }

                            // Get photons
                            GPhotons photons = model->mc(m_area, dir, rad,
                                                         emin, emax,
                                                         tstart, tstop, m_ran);

                            // Simulate events from photons
                            for (int i = 0; i < photons.size(); ++i) {

                                // Apply deadtime correction
                                if (m_deadc < 1.0) {
                                    if (m_ran.uniform() > m_deadc) {
                                        ndeadc++;
                                        continue;
                                    }
                                }

                                // Increment photon counter
                                nphotons++;

                                // Simulate event
                                GCTAEventAtom* event = rsp->mc(m_area, photons[i], *pnt, m_ran);
                                if (event != NULL) {
                                    events->append(*event);
                                    delete event;
                                }

                            } // endfor: looped over events

                            // Go to next time slice
                            tstart = tstop;
                            tstop  = tstart + m_time_max;
            
                            // Reset indentation
                            if (logExplicit()) {
                                if (tmax - tmin > m_time_max) {
                                    indent--;
                                    log.indent(indent);
                                }
                            }
            
                        } // endwhile: looped over time slices

                        // Dump simulation results
                        if (logNormal()) {
                            log << parformat("MC source photons");
                            log << str(nphotons);
                            if (model->name().length() > 0) {
                                log << " [" << model->name() << "]";
                            }
                            log << std::endl;
                            log << parformat("MC photons during deadtime");
                            log << str(ndeadc);
                            if (model->name().length() > 0) {
                                log << " [" << model->name() << "]";
                            }
                            log << std::endl;
                            log << parformat("MC source events");
                            log << str(events->size());
                            if (model->name().length() > 0) {
                                log << " [" << model->name() << "]";
                            }
                            log << std::endl;
                        }

                    } // endif: model was a sky model

                } // endfor: looped over models

                // Dump simulation results
                if (logNormal()) {
                    log << parformat("MC source events");
                    log << str(events->size());
                    log << " (all source models)";
                    log << std::endl;
                }
                    
                // Reset indentation
                if (logNormal()) {
                    if (events->ebounds().size() > 1) {
                        indent--;
                        log.indent(indent);
                    }
                }
                
            } // endfor: looped over all energy boundaries

            // Reset indentation
            if (logNormal()) {
                if (events->gti().size() > 1) {
                    indent--;
                    log.indent(indent);
                }
            }

        } // endfor: looped over all time intervals

        // Reset indentation
        log.indent(0);

    } // endif: observation pointer was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Simulate background events from model
 *
 * @param[in] obs Pointer on CTA observation.
 * @param[in] models Models.
 *
 * Simulate background events from models. The events are stored as event
 * list in the observation.
 *
 * This method does nothing if the observation pointer is NULL.
 ***************************************************************************/
void ctobssim::simulate_background(GCTAObservation* obs, const GModels& models)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Make sure that the observation holds a CTA event list. If this
        // is not the case then allocate and attach a CTA event list now.
        if (dynamic_cast<const GCTAEventList*>(obs->events()) == NULL) {
            set_list(obs);
        }

        // Get pointer on event list (circumvent const correctness)
        GCTAEventList* events = (GCTAEventList*)(obs->events());

        // Loop over all models
        for (int i = 0; i < models.size(); ++i) {

            // Get model (NULL if not a radial acceptance model)
            const GCTAModelRadialAcceptance* model = 
                  dynamic_cast<const GCTAModelRadialAcceptance*>(&models[i]);

            // If we have a radial acceptance model then simulate events
            if (model != NULL) {

                // Get simulated event list
                GCTAEventList* list = model->mc(*obs, m_ran);

                // Reserves space for events
                events->reserve(list->size()+events->size());

                // Initialise event counters
                int nevents = 0;
                int ndeadc  = 0;

                // Append events
                for (int k = 0; k < list->size(); k++) {

                    // Apply deadtime correction
                    if (m_deadc < 1.0) {
                        if (m_ran.uniform() > m_deadc) {
                            ndeadc++;
                            continue;
                        }
                    }

                    // Increment background event counter
                    nevents++;

                    // Append event
                    events->append(*((*list)[k]));

                } // endfor: looped over all events

                // Dump simulation results
                if (logNormal()) {
                    log << parformat("MC background events");
                    log << str(nevents) << std::endl;
                    log << parformat("MC events during deadtime");
                    log << str(ndeadc) << std::endl;
                }

                // Free event list
                delete list;

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
    m_obs.clear();

    // Set fixed parameters
    m_area = 19634954.0 * 1.0e4;     //!< pi*(2500^2) m^2
    m_time_max.met(1800.0);          //!< Maximum length of time slice (s)

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
    m_area     = app.m_area;
    m_time_max = app.m_time_max;
    m_ran      = app.m_ran;
    m_obs      = app.m_obs;

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
