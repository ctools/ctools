/***************************************************************************
 *                        ctbin - Event binning tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
#define G_INIT_CUBE                                      "ctbin::init_cube()"
#define G_FILL_CUBE                      "ctbin::fill_cube(GCTAObservation*)"
#define G_GET_EBOUNDS                                  "ctbin::get_ebounds()"

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
 * This method creates an instance of the class by copying an existing
 * observations container.
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
 * @param[in] app Application.
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
 * @brief Clear instance
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
 * @brief Execute application
 *
 * This is the main execution method of the ctbin class. It is invoked when
 * the executable is called from command line.
 *
 * The method reads the task parameters, bins the event list(s) into counts
 * map(s), and writes the results into FITS files on disk.
 ***************************************************************************/
void ctbin::execute(void)
{
    // Signal that some parameters should be read ahead
    m_read_ahead = true;

    // Bin the event data
    run();

    // Save the counts map into FITS file
    save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Bin the event data
 *
 * This method loops over all observations found in the observation conatiner
 * and bins all events from the event list(s) into counts map(s). Note that
 * each event list is binned in a separate counts map, hence no summing of
 * events is done.
 ***************************************************************************/
void ctbin::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Write parameters into logger
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
            log.header1("Bin observations");
        }
        else {
            log.header1("Bin observation");
        }
    }

    // Initialise counts cube
    init_cube();

    // Loop over all observations in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

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

            // Fill the cube
            fill_cube(obs);

        } // endif: CTA observation found

    } // endfor: looped over observations

    // Set a single cube in the observation container
    obs_cube();

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Binned observations");
        }
        else {
            log.header1("Binned observation");
        }
        log << m_obs << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save counts cube
 *
 * This method saves the counts cube into a FITS file.
 ***************************************************************************/
void ctbin::save(void)
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

    // Get CTA observation from observation container
    GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

    // Save only if observation is valid
    if (obs != NULL) {
        obs->save(m_outfile, clobber());
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
 ***************************************************************************/
void ctbin::init_members(void)
{
    // Initialise members
    m_evfile.clear();
    m_outfile.clear();
    m_proj.clear();
    m_coordsys.clear();
    m_usepnt   = false;
    m_xref     = 0.0;
    m_yref     = 0.0;
    m_binsz    = 0.0;
    m_nxpix    = 0;
    m_nypix    = 0;

    // Initialise protected members
    m_obs.clear();
    m_read_ahead = false;
    m_cube.clear();
    m_ebounds.clear();
    m_gti.clear();
    m_ontime   = 0.0;
    m_livetime = 0.0;

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
void ctbin::copy_members(const ctbin& app)
{
    // Copy attributes
    m_evfile   = app.m_evfile;
    m_outfile  = app.m_outfile;
    m_proj     = app.m_proj;
    m_coordsys = app.m_coordsys;
    m_usepnt   = app.m_usepnt;
    m_xref     = app.m_xref;
    m_yref     = app.m_yref;
    m_binsz    = app.m_binsz;
    m_nxpix    = app.m_nxpix;
    m_nypix    = app.m_nypix;

    // Copy protected members
    m_obs        = app.m_obs;
    m_read_ahead = app.m_read_ahead;
    m_cube       = app.m_cube;
    m_ebounds    = app.m_ebounds;
    m_gti        = app.m_gti;
    m_ontime     = app.m_ontime;
    m_livetime   = app.m_livetime;

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
    // Initialize energy boundaries
    m_ebounds.clear();

    // If there are no observations in container then add a single CTA
    // observation using the parameters from the parameter file
    if (m_obs.size() == 0) {

        // Get name of CTA events file
        m_evfile = (*this)["evfile"].filename();

        // Allocate CTA observation
        GCTAObservation obs;

        // Try first to open as FITS file
        try {

            // Load event list in CTA observation
            obs.load(m_evfile);

            // Append CTA observation to container
            m_obs.append(obs);

        }

        // ... otherwise try to open as XML file
        catch (GException::fits_open_error &e) {

            // Load observations from XML file
            m_obs.load(m_evfile);

        }

    } // endif: there was no observation in the container

    // Get parameters
    m_usepnt = (*this)["usepnt"].boolean();
    if (!m_usepnt) {
        m_xref = (*this)["xref"].real();
        m_yref = (*this)["yref"].real();
    }
    m_proj     = (*this)["proj"].string();
    m_coordsys = (*this)["coordsys"].string();
    m_binsz    = (*this)["binsz"].real();
    m_nxpix    = (*this)["nxpix"].integer();
    m_nypix    = (*this)["nypix"].integer();
    
    // Set energy boundaries
    m_ebounds  = get_ebounds();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (m_read_ahead) {
        m_outfile = (*this)["outfile"].filename();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise counts cube information
 *
 * @exception GException::invalid_value
 *            No valid CTA observation found to derive the counts cube
 *            map centre.
 *
 * Initialises the skymap, energy boundaries and GTI for a counts cube.
 ***************************************************************************/
void ctbin::init_cube(void)
{
    // Initialse cube information
    m_ontime   = 0.0;
    m_livetime = 0.0;
    m_cube.clear();
    m_gti.clear();

    // Set event cube centre, either from the user parameters or from the
    // pointing
    double xref = m_xref;
    double yref = m_yref;
    if (m_usepnt) {

        // Dummy: get pointing from first observation. Ultimately, we want
        // to get the pointing from a kind of average
        bool found = false;
        for (int i = 0; i < m_obs.size(); ++i) {
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
            if (obs != NULL) {
                const GCTAPointing& pnt = obs->pointing();
                if (gammalib::toupper(m_coordsys) == "GAL") {
                    xref = pnt.dir().l_deg();
                    yref = pnt.dir().b_deg();
                }
                else {
                    xref = pnt.dir().ra_deg();
                    yref = pnt.dir().dec_deg();
                }
                found = true;
                break;
            }
        }

        // Signal if no pointing is found
        if (!found) {
            std::string msg = "No valid CTA observation has been found in "
                              "observation list, hence no pointing information "
                              "could be extracted. Use the \"usepnt=no\" "
                              "option and specify pointing explicitly when "
                              "running ctbin.";
            throw GException::invalid_value(G_INIT_CUBE, msg);
        }

    } // endif: used pointing

    // Create skymap
    m_cube = GSkymap(m_proj, m_coordsys,
                     xref, yref, -m_binsz, m_binsz,
                     m_nxpix, m_nypix, m_ebounds.size());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill events into counts cube
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::invalid_value
 *            No event list found in observation.
 *
 * Fills the events from an event list in the counts cube setup by init_cube.
 ***************************************************************************/
void ctbin::fill_cube(GCTAObservation* obs)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Make sure that the observation holds a CTA event list. If this
        // is not the case then throw an exception.
        const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(obs->events());
        if (events == NULL) {
            std::string msg = "CTA Observation does not contain an event "
                              "list. Event list information is needed to "
                              "fill the counts map.";
            throw GException::invalid_value(G_FILL_CUBE, msg);
        }

        // Get the RoI
        const GCTARoi& roi = events->roi();

        // Initialise binning statistics
        int num_outside_roi  = 0;
        int num_outside_map  = 0;
        int num_outside_ebds = 0;
        int num_in_map       = 0;

        // Fill sky map
        for (int i = 0; i < events->size(); ++i) {

            // Get event
            const GCTAEventAtom* event = (*events)[i];

            // Determine sky pixel
            GCTAInstDir* inst  = (GCTAInstDir*)&(event->dir());
            GSkyDir      dir   = inst->dir();
            GSkyPixel    pixel = m_cube.dir2pix(dir);

            // Skip if pixel is outside RoI
            if (roi.centre().dir().dist_deg(dir) > roi.radius()) {
                num_outside_roi++;
                continue;
            }

            // Skip if pixel is out of range
            if (pixel.x() < -0.5 || pixel.x() > (m_nxpix-0.5) ||
                pixel.y() < -0.5 || pixel.y() > (m_nypix-0.5)) {
                num_outside_map++;
                continue;
            }

            // Determine energy bin. Skip if we are outside the energy range
            int index = m_ebounds.index(event->energy());
            if (index == -1) {
                num_outside_ebds++;
                continue;
            }

            // Fill event in skymap
            m_cube(pixel, index) += 1.0;
            num_in_map++;

        } // endfor: looped over all events

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
            log << m_cube << std::endl;
        }

    } // endif: observation was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create output observation container.
 *
 * Creates an output observation container that combines all input CTA
 * observation into a single cube-style observation. All non CTA observations
 * present in the observation container are kept. The method furthermore
 * conserves any response information in case that a single CTA observation
 * is provided. This supports the original binned analysis.
 ***************************************************************************/
void ctbin::obs_cube(void)
{
    // If we have only a single CTA observation in the container, then
    // keep that observation and just attach the event cube to it
    if (m_obs.size() == 1) {

        // Attach event cube to CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);
        if (obs != NULL) {
            obs->events(this->cube());
        }

    }

    // ... otherwise put a single CTA observation in container
    else {

        // Allocate observation container
        GObservations container;

        // Allocate CTA observation.
        GCTAObservation obs;

        // Attach event cube to CTA observation
        obs.events(this->cube());

        // Set map centre as pointing
        GSkyPixel    pixel(0.5*double(m_cube.nx()), 0.5*double(m_cube.ny()));
        GSkyDir      centre = m_cube.pix2dir(pixel);
        GCTAPointing pointing(centre);

        // Compute deadtime correction
        double deadc = (m_ontime > 0.0) ? m_livetime / m_ontime : 0.0;

        // Set CTA observation attributes
        obs.pointing(pointing);
        obs.obs_id(0);
        obs.ra_obj(centre.ra_deg());   //!< Dummy
        obs.dec_obj(centre.dec_deg()); //!< Dummy
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
        }

        // Set observation container
        m_obs = container;

    } // endelse: there was not a single CTA observation

    // Return
    return;
}
