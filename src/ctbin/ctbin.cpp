/***************************************************************************
 *                      ctbin - CTA data binning tool                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @brief CTA data binning tool implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctbin.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_BIN_EVENTS                    "ctbin::bin_events(GCTAObservation*)"

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
ctbin::ctbin(void) : GApplication(CTBIN_NAME, CTBIN_VERSION)
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
 * This method creates an instance of the class by copying an existing
 * observations container.
 ***************************************************************************/
ctbin::ctbin(GObservations obs) : GApplication(CTBIN_NAME, CTBIN_VERSION)
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
ctbin::ctbin(int argc, char *argv[]) : 
                          GApplication(CTBIN_NAME, CTBIN_VERSION, argc, argv)
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
ctbin::ctbin(const ctbin& app) : GApplication(app)
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
ctbin& ctbin::operator= (const ctbin& app)
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
void ctbin::clear(void)
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
 * This method bins the events data into a counts maps and saves the counts
 * map into a FITS file.
 ***************************************************************************/
void ctbin::execute(void)
{
    // Read ahead output filename so that it gets dumped correctly in the
    // parameters log
    m_outfile = (*this)["outfile"].filename();

    // Bin the event data
    run();

    // Save the counts map into FITS file
    save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Bin the event data
 ***************************************************************************/
void ctbin::run(void)
{
    // Switch screen logging on in debug mode
    if (logDebug())
        log.cout(true);

    // Get parameters
    get_parameters();

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1)
            log.header1("Observations");
        else
            log.header1("Observation");
        log << m_obs << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1)
            log.header1("Bin observations");
        else
            log.header1("Bin observation");
    }

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(&m_obs[i]);

        // Continue only if observation is a CTA observation
        if (obs != NULL) {

            // Write header for observation
            if (logTerse()) {
                if (obs->name().length() > 1)
                    log.header3("Observation "+obs->name());
                else
                    log.header3("Observation");
            }

            // Bin events into counts map
            bin_events(obs);

        } // endif: CTA observation found

    } // endfor: looped over observations

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1)
            log.header1("Binned observations");
        else
            log.header1("Binned observation");
        log << m_obs << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save observation
 *
 * This method saves the counts map(s) into (a) FITS file(s).
 ***************************************************************************/
void ctbin::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1)
            log.header1("Save observations");
        else
            log.header1("Save observation");
    }

    // Get output filename
    m_outfile = (*this)["outfile"].filename();

    // Loop over all observations in the container
    int file_num = 0;
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(&m_obs[i]);

        // Save only if observation is a CTA observation
        if (obs != NULL) {

            // Set filename. If more than one file will be created an
            // index "_xxx" will be appended.
            std::string filename = m_outfile;
            if (file_num > 0)
                filename += "_"+str(file_num);

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
void ctbin::get_parameters(void)
{
    // If there are no observations in container then add a single CTA
    // observation using the parameters from the parameter file
    if (m_obs.size() == 0) {

        // Get name of CTA events file
        m_evfile = (*this)["evfile"].filename();

        // Load unbinned CTA observation
        GCTAObservation obs;
        obs.load_unbinned(m_evfile);

        // Append CTA observation to container
        m_obs.append(obs);

    } // endif: there was no observation in the container

    // Get remaining parameters
    m_emin     = (*this)["emin"].real();
    m_emax     = (*this)["emax"].real();
    m_enumbins = (*this)["enumbins"].integer();
    m_proj     = (*this)["proj"].string();
    m_coordsys = (*this)["coordsys"].string();
    m_xref     = (*this)["xref"].real();
    m_yref     = (*this)["yref"].real();
    m_binsz    = (*this)["binsz"].real();
    m_nxpix    = (*this)["nxpix"].integer();
    m_nypix    = (*this)["nypix"].integer();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Bin events into a counts map
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::no_list
 *            No event list found in observation.
 *
 * This method bins the events found in a CTA events list into a counts map
 * and replaces the event list by the counts map in the observation. The
 * energy boundaries of the counts map are also stored in the observation's
 * energy boundary member.
 ***************************************************************************/
void ctbin::bin_events(GCTAObservation* obs)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Make sure that the observation holds a CTA event list. If this
        // is not the case then throw an exception.
        if (dynamic_cast<const GCTAEventList*>(obs->events()) == NULL)
            throw GException::no_list(G_BIN_EVENTS);

        // Setup energy range covered by data
        GEnergy  emin;
        GEnergy  emax;
        GEbounds ebds;
        emin.TeV(m_emin);
        emax.TeV(m_emax);
        ebds.setlog(emin, emax, m_enumbins);

        // Get Good Time intervals
        GGti gti = obs->events()->gti();

        // Create skymap
        GSkymap map = GSkymap(m_proj, m_coordsys,
                             m_xref, m_yref, m_binsz, m_binsz,
                             m_nxpix, m_nypix, m_enumbins);

        // Initialise binning statistics
        int num_outside_map  = 0;
        int num_outside_ebds = 0;
        int num_in_map       = 0;

        // Fill sky map
        GCTAEventList* events = (GCTAEventList*)obs->events();
        for (GCTAEventList::iterator event = events->begin(); event != events->end(); ++event) {

            // Determine sky pixel
            GCTAInstDir* inst  = (GCTAInstDir*)&(event->dir());
            GSkyDir      dir   = inst->skydir();
            GSkyPixel    pixel = map.dir2xy(dir);

            // Skip if pixel is out of range
            if (pixel.x() < -0.5 || pixel.x() > (m_nxpix-0.5) ||
                pixel.y() < -0.5 || pixel.y() > (m_nypix-0.5)) {
                num_outside_map++;
                continue;
            }

            // Determine energy bin. Skip if we are outside the energy range
            int index = ebds.index(event->energy());
            if (index == -1) {
                num_outside_ebds++;
                continue;
            }

            // Fill event in skymap
            map(pixel, index) += 1.0;
            num_in_map++;

        } // endfor: looped over all events

        // Log binning results
        if (logTerse()) {
            log << std::endl;
            log.header1("Binning");
            log << parformat("Events in list");
            log << obs->events()->size() << std::endl;
            log << parformat("Events in map");
            log << num_in_map << std::endl;
            log << parformat("Events outside map area");
            log << num_outside_map << std::endl;
            log << parformat("Events outside energy bins");
            log << num_outside_ebds << std::endl;
        }

        // Log map
        if (logTerse()) {
            log << std::endl;
            log.header1("Counts map");
            log << map << std::endl;
        }

        // Create events cube from sky map
        GCTAEventCube cube(map, ebds, gti);

        // Replace event list by event cube in observation
        obs->events(&cube);

    } // endif: observation was valid

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
    m_obs.clear();
    m_emin     = 0.0;
    m_emax     = 0.0;
    m_enumbins = 0;
    m_xref     = 0.0;
    m_yref     = 0.0;
    m_binsz    = 0.0;
    m_nxpix    = 0;
    m_nypix    = 0;

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
    m_obs      = app.m_obs;
    m_emin     = app.m_emin;
    m_emax     = app.m_emax;
    m_enumbins = app.m_enumbins;
    m_xref     = app.m_xref;
    m_yref     = app.m_yref;
    m_binsz    = app.m_binsz;
    m_nxpix    = app.m_nxpix;
    m_nypix    = app.m_nypix;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctbin::free_members(void)
{
    // Write separator into logger
    if (logTerse())
        log << std::endl;

    // Return
    return;
}
