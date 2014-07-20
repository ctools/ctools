/***************************************************************************
 *                     ctskymap - CTA sky mapping tool                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @brief CTA sky mapping tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctskymap.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_INIT_MAP                 "ctskymap::init_map(GCTAObservation* obs)"
#define G_BIN_EVENTS                 "ctskymap::bin_events(GCTAObservation*)"

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
ctskymap::ctskymap(void) : GApplication(CTSKYMAP_NAME, CTSKYMAP_VERSION)
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
 * param[in] obs Observation container.
 *
 * This method creates an instance of the class by copying an existing
 * observations container.
 ***************************************************************************/
ctskymap::ctskymap(const GObservations& obs) :
          GApplication(CTSKYMAP_NAME, CTSKYMAP_VERSION)
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
ctskymap::ctskymap(int argc, char *argv[]) : 
                    GApplication(CTSKYMAP_NAME, CTSKYMAP_VERSION, argc, argv)
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
ctskymap::ctskymap(const ctskymap& app) : GApplication(app)
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
 * @param[in] app Application.
 ***************************************************************************/
ctskymap& ctskymap::operator= (const ctskymap& app)
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
void ctskymap::clear(void)
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
 * This method creates the sky map and saves the map into a FITS file.
 ***************************************************************************/
void ctskymap::execute(void)
{
    // Read ahead output filename so that it gets dumped correctly in the
    // parameters log
    m_outfile = (*this)["outfile"].filename();

    // Create the sky map
    run();

    // Save the sky map into FITS file
    save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Creates sky maps from data
 *
 * This method is the main code. It
 * (1) reads task parameters from the par file
 * (2) initialises the sky maps
 * (3) loops over all observations to add its events to the sky map
 ***************************************************************************/
void ctskymap::run(void)
{
    // Initialise statistics
    int num_obs = 0;

    // Switch screen logging on in debug mode
    if (logDebug()) {
        log.cout(true);
    }

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
            log.header1("Map observations");
        }
        else {
            log.header1("Map observation");
        }
    }

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Continue only if observation is a CTA observation
        if (obs != NULL) {

            // Write header for current observation
            if (logTerse()) {
                if (obs->name().length() > 1) {
                    log.header3("Observation "+obs->name());
                }
                else {
                    log.header3("Observation");
                }
            }

            // If this is the first valid observation we're working on
            // then initialise the sky maps
            if (num_obs == 0) {
                init_map(obs);
            }

            // Map events into sky map
            map_events(obs);
            
            // Increment observation counter
            num_obs++;

        } // endif: CTA observation found

    } // endfor: looped over observations

    // Write observation(s) into logger
/*
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Map observations");
        }
        else {
            log.header1("Map observation");
        }
        log << m_obs << std::endl;
    }
*/

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save observation
 *
 * This method saves the counts map(s) into (a) FITS file(s).
 ***************************************************************************/
void ctskymap::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save sky map");
    }

    // Get output filename
    m_outfile = (*this)["outfile"].filename();

    // Save sky map
    m_skymap.save(m_outfile, clobber());

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
    // If there are no observations in container then add a single CTA
    // observation using the parameters from the parameter file
    if (m_obs.size() == 0) {

        // Get name of CTA events file
        m_evfile = (*this)["evfile"].filename();

        // Load unbinned CTA observation
        GCTAObservation obs;
        obs.load(m_evfile);

        // Append CTA observation to container
        m_obs.append(obs);
        
        // Use the xref and yref parameters for binning (otherwise the
        // pointing direction(s) is/are used)
        m_xref = (*this)["xref"].real();
        m_yref = (*this)["yref"].real();

    } // endif: there was no observation in the container

    // Get remaining parameters
    m_emin     = (*this)["emin"].real();
    m_emax     = (*this)["emax"].real();
    m_proj     = (*this)["proj"].string();
    m_coordsys = (*this)["coordsys"].string();
    m_binsz    = (*this)["binsz"].real();
    m_nxpix    = (*this)["nxpix"].integer();
    m_nypix    = (*this)["nypix"].integer();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise sky map
 *
 * @param[in] obs CTA observation (no NULL pointer allowed).
 *
 * This method initialises the sky map.
 ***************************************************************************/
void ctskymap::init_map(GCTAObservation* obs)
{
    // Clear any existing sky map
    m_skymap.clear();
    
    // Get map centre. If no centre is specified then extract the map centre
    // from the pointing direction. This obviously only works for 
    double xref;
    double yref;
    if (m_xref != 9999.0 && m_yref != 9999.0) {
        xref = m_xref;
        yref = m_yref;
    }
    else {


        // Get pointer on CTA pointing
        const GCTAPointing& pnt = obs->pointing();
            
        // Set reference point to pointing
        if (gammalib::toupper(m_coordsys) == "GAL") {
            xref = pnt.dir().l_deg();
            yref = pnt.dir().b_deg();
        }
        else {
            xref = pnt.dir().ra_deg();
            yref = pnt.dir().dec_deg();
        }

    } // endelse: map centre set to pointing

    // Create skymap
    m_skymap = GSkymap(m_proj, m_coordsys, xref, yref, -m_binsz, m_binsz,
                       m_nxpix, m_nypix, 1);

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
 * @exception GCTAException::no_pointing
 *            No valid CTA pointing found.
 *
 * This method maps the events found in a CTA events list into a sky map.
 ***************************************************************************/
void ctskymap::map_events(GCTAObservation* obs)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Make sure that the observation holds a CTA event list. If this
        // is not the case then throw an exception.
        if (dynamic_cast<const GCTAEventList*>(obs->events()) == NULL) {
            throw GException::no_list(G_BIN_EVENTS);
        }

        // Setup energy range covered by data
        GEnergy  emin;
        GEnergy  emax;
        GEbounds ebds;
        emin.TeV(m_emin);
        emax.TeV(m_emax);

        // Initialise binning statistics
        int num_outside_map    = 0;
        int num_outside_erange = 0;
        int num_in_map         = 0;

        // Fill sky map
        GCTAEventList* events = static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));
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
            if (pixel.x() < -0.5 || pixel.x() > (m_nxpix-0.5) ||
                pixel.y() < -0.5 || pixel.y() > (m_nypix-0.5)) {
                num_outside_map++;
                continue;
            }

            // Fill event in skymap
            m_skymap(pixel, 0) += 1.0;
            num_in_map++;

        } // endfor: looped over all events

        // Log binning results
        if (logTerse()) {
            log << std::endl;
            log.header1("Mapping");
            log << gammalib::parformat("Events in list");
            log << obs->events()->size() << std::endl;
            log << gammalib::parformat("Events in map");
            log << num_in_map << std::endl;
            log << gammalib::parformat("Events outside map area");
            log << num_outside_map << std::endl;
            log << gammalib::parformat("Events outside energy range");
            log << num_outside_erange << std::endl;
        }

        // Log map
        if (logTerse()) {
            log << std::endl;
            log.header1("Sky map");
            log << m_skymap << std::endl;
        }

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
void ctskymap::init_members(void)
{
    // Initialise members
    m_evfile.clear();
    m_outfile.clear();
    m_proj.clear();
    m_coordsys.clear();
    m_obs.clear();
    m_skymap.clear();
    m_emin     = 0.0;
    m_emax     = 0.0;
    m_xref     = 9999.0; // Flags unset (use pointing direction)
    m_yref     = 9999.0; // Flags unset (use pointing direction)
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
void ctskymap::copy_members(const ctskymap& app)
{
    // Copy attributes
    m_evfile   = app.m_evfile;
    m_outfile  = app.m_outfile;
    m_proj     = app.m_proj;
    m_coordsys = app.m_coordsys;
    m_obs      = app.m_obs;
    m_skymap   = app.m_skymap;
    m_emin     = app.m_emin;
    m_emax     = app.m_emax;
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
void ctskymap::free_members(void)
{
    // Write separator into logger
    if (logTerse())
        log << std::endl;

    // Return
    return;
}
