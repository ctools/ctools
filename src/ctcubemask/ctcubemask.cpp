/***************************************************************************
 *                      ctcubemask - Cube filter tool                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Chia-Chun Lu                                *
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
 * @file ctcubemask.cpp
 * @brief Cube filter tool implementation
 * @author Chia-Chun Lu
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctcubemask.hpp"
#include "GTools.hpp"


/* __ Method name definitions ____________________________________________ */
#define G_RUN                                             "ctcubemask::run()"
#define G_APPLY_MASK               "ctcubemask::apply_mask(GCTAObservation*)"

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
ctcubemask::ctcubemask(void) : ctobservation(CTCUBEMASK_NAME, CTCUBEMASK_VERSION)
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
 * This constructor creates an instance of the class that is initialised from
 * an observation container.
 ***************************************************************************/
ctcubemask::ctcubemask(const GObservations& obs) :
            ctobservation(CTCUBEMASK_NAME, CTCUBEMASK_VERSION, obs)
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
 ***************************************************************************/
ctcubemask::ctcubemask(int argc, char *argv[]) : 
            ctobservation(CTCUBEMASK_NAME, CTCUBEMASK_VERSION, argc, argv)
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
ctcubemask::ctcubemask(const ctcubemask& app) : ctobservation(app)
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
ctcubemask::~ctcubemask(void)
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
ctcubemask& ctcubemask::operator=(const ctcubemask& app)
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
 * @brief Clear ctcubemask tool
 *
 * Clears ctcubemask tool.
 ***************************************************************************/
void ctcubemask::clear(void)
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
 * @brief mask data cube
 *
 * This method reads in the application parameters and loops over all
 * observations that were found to apply a mask on the event cube.
 ***************************************************************************/
void ctcubemask::run(void)
{
    // Switch screen logging on in debug mode
    if (logDebug()) {
        log.cout(true);
    }

    // Get parameters
    get_parameters();

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

    // Write header into logger
    log_header1(TERSE, "Apply mask");

    // Initialise counters
    int n_observations = 0;

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Initialise event input and output filenames
        m_infiles.push_back("");

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Continue only if observation is a CTA observation
        if (obs != NULL) {

            // Write header for the current observation
            log_header3(TERSE, get_obs_header(obs));

            // Increment counter
            n_observations++;

            // Save event file name (for possible saving)
            m_infiles[i] = obs->eventfile();

            // Apply mask on the event cube
            apply_mask(obs);

        } // endif: had a CTA observation

    } // endfor: looped over all observations

    // If more than a single observation has been handled then make sure that
    // an XML file will be used for storage
    if (n_observations > 1) {
        m_use_xml = true;
    }

    // Write resulting observation container into logger
    log_observations(NORMAL, m_obs, "Observations after event bin masking");

    // Optionally publish counts cube
    if (m_publish) {
        publish();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save the masked event cube(s)
 *
 * This method saves the masked event cube(s) into FITS file(s). There are
 * two modes, depending on the m_use_xml flag.
 *
 * If m_use_xml is true, all masked event cube(s) will be saved into FITS
 * files, where the output filenames are constructued from the input
 * filenames by prepending the m_prefix string to name. Any path information
 * will be stripped form the input name, hence event cube files will be written
 * into the local working directory (unless some path information is present
 * in the prefix). In addition, an XML file will be created that gathers
 * the filename information for the masked event cube(s). If an XML file
 * was present on input, all metadata information will be copied from this
 * input file.
 *
 * If m_use_xml is false, the masked event cubes will be saved into a FITS
 * file.
 ***************************************************************************/
void ctcubemask::save(void)
{
    // Write header into logger
    log_header1(TERSE, gammalib::number("Save counts cube", m_obs.size()));

    // Get counts cube filename
    m_outcube = (*this)["outcube"].filename();

    // Case A: Save counts cube(s) and XML metadata information
    if (m_use_xml) {

        // Get prefix
        m_prefix = (*this)["prefix"].string();

        // Save XML
        save_xml();
    }

    // Case B: Save counts cube as FITS file
    else {
        save_fits();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Publish counts cube
 *
 * @param[in] name Counts cube name.
 *
 * Publishes the first counts cube in the observation container on the VO
 * Hub.
 ***************************************************************************/
void ctcubemask::publish(const std::string& name)
{
    // Write header into logger
    log_header1(TERSE, "Publish counts cube");

    // Set default name is user name is empty
    std::string user_name(name);
    if (user_name.empty()) {
        user_name = CTCUBEMASK_NAME;
    }

    // Get first CTA observation from observation container
    GCTAObservation* obs  = dynamic_cast<GCTAObservation*>(m_obs[0]);
    if (obs != NULL) {

        // Get counts cube
        GCTAEventCube* cube = dynamic_cast<GCTAEventCube*>(obs->events());
        if (cube != NULL) {

            // Log counts cube name
            log_value(NORMAL, "Counts cube name", user_name);

            // Publish counts cube
            cube->counts().publish(user_name);

        }
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
void ctcubemask::init_members(void)
{
    // Initialise parameters
	m_regfile.clear();
    m_outcube.clear();
	m_prefix.clear();
    m_usepnt  = false;
    m_roi.clear();
    m_emin    = 0.0;
    m_emax    = 0.0;
    m_publish = false;

    // Initialise protected members
    m_infiles.clear();
    m_select_energy = true;
    m_select_roi    = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctcubemask::copy_members(const ctcubemask& app)
{
    // Copy parameters
    m_regfile = app.m_regfile;
    m_outcube = app.m_outcube;
    m_prefix  = app.m_prefix;
    m_usepnt  = app.m_usepnt;
    m_roi     = app.m_roi;
    m_emin    = app.m_emin;
    m_emax    = app.m_emax;
    m_publish = app.m_publish;

    // Copy protected members
    m_infiles       = app.m_infiles;
    m_select_energy = app.m_select_energy;
    m_select_roi    = app.m_select_roi;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctcubemask::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. 
 *
 * This method also loads observations if no observations are yet allocated.
 * Observations are either loaded from a single CTA even cube, or from a
 * XML file using the metadata information that is stored in that file.
 ***************************************************************************/
void ctcubemask::get_parameters(void)
{
    // Initialise selection flags
    m_select_energy = true;
    m_select_roi    = true;

    // Setup observations from "inobs" parameter. Do not request response
    // information and do not accept event lists.
    setup_observations(m_obs, false, false, true);

    // Get parameters
    m_regfile = (*this)["regfile"].filename();
    m_usepnt  = (*this)["usepnt"].boolean();

    // Get the RoI and enable RoI selection if the RoI is valid
    m_roi        = get_roi();
    m_select_roi = m_roi.is_valid();

    // Check for sanity of energy selection parameters
    if ((*this)["emin"].is_valid() && (*this)["emax"].is_valid()) {
        m_emin          = (*this)["emin"].real();
        m_emax          = (*this)["emax"].real();
        m_select_energy = true;
    }
    else {
        m_select_energy = false;
    }

    // Get remaining parameters
    m_publish = (*this)["publish"].boolean();

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (read_ahead()) {
        m_outcube = (*this)["outcube"].filename();
        m_prefix  = (*this)["prefix"].string();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Apply mask to event cube
 *
 * @param[in] obs CTA observation.
 *
 * Apply mask to an event cube. The mask sets content of event bins 
 * outside the region of interest or the energy of interest to -1. 
 * These pixels will be ignored in likelihood fitting.
 ***************************************************************************/
void ctcubemask::apply_mask(GCTAObservation* obs)
{
    // Get pointer to CTA event cube. Continue only if we have a cube
    GCTAEventCube* cube = dynamic_cast<GCTAEventCube*>
                                      (const_cast<GEvents*>(obs->events()));
    if (cube != NULL) {

        // Extract event cube and energy boundaries
        GSkyMap         map     = cube->counts();
        const GEbounds& ebounds = cube->ebounds();

        // If no energy selection is required set energy boundaries to cube boundaries
        if (!m_select_energy) {
            m_emin = ebounds.emin().TeV();
            m_emax = ebounds.emax().TeV();
        }

        // Initialise energy selection
        int npix   = map.npix();
        int n_ebin = ebounds.size();
        int e_idx1 = 0;
        int e_idx2 = n_ebin;

        // Determine number of events before masking
        double sum_before = 0.0;
        for (int i = 0; i < n_ebin; ++i) {
            for (int pixel = 0; pixel < npix; ++pixel) {
                if (map(pixel, i) >= 0.0) {
                    sum_before += map(pixel, i);
                }
            }
        }

        // Loop over ebounds to find the first valid energy band
        for (int i = 0; i < n_ebin; ++i) {
            double emin = ebounds.emin(i).TeV() + 1.0e-6; // Rounding tolerance
            if (emin >= m_emin) {
                e_idx1 = i;
                break;
            }
        }

        // Loop over ebounds to find the last valid energy band
        for (int i = n_ebin-1; i >= 0; i--) {
            double emax = ebounds.emax(i).TeV() - 1.0e-6; // Rounding tolerance
            if (emax <= m_emax) {
                e_idx2 = i;
                break;
            }
        }

        // Set all pixels outside the desired energy bands to -1.0
        for (int i = 0; i < e_idx1; ++i) {
            for (int pixel = 0; pixel < npix; ++pixel) {
                map(pixel,i) = -1.0;
            }
        }
        for (int i = e_idx2 + 1; i < n_ebin; ++i) {
            for (int pixel = 0; pixel < npix; ++pixel) {
                map(pixel,i) = -1.0;
            }
        }

        // Log selected energy band
        log_value(NORMAL, "Selected energy band",
                          gammalib::str(ebounds.emin(e_idx1).TeV())+" - "+
                          gammalib::str(ebounds.emax(e_idx2).TeV())+" TeV");

        // Set all pixels inside selected energy bands but outside RoI
        // to -1.0 if requested
        if (m_select_roi) {
            GSkyRegionCircle roi(m_roi.centre().dir(), m_roi.radius());
            for (int i = e_idx1; i <= e_idx2; ++i) {
                for (int pixel = 0; pixel < npix; ++pixel) {
                    GSkyDir dir = map.inx2dir(pixel);
                    if (!roi.contains(dir)) {
                        map(pixel,i) = -1.0;
                    }
                }
            }

            // Log selected RoI
            log_value(NORMAL, "Selected RoI",
                      "RA="+gammalib::str(roi.centre().ra_deg())+" deg, "+
                      "DEC="+gammalib::str(roi.centre().dec_deg())+" deg, "+
                      "Radius="+gammalib::str(roi.radius())+" deg");

        } // endif: applied RoI selection

        // Set all pixels inside selected energy bands and inside exclusion
        // regions to -1.0
        if (m_regfile != "NONE") {
            GSkyRegions regions = GSkyRegions(m_regfile);
            for (int i = e_idx1; i <= e_idx2; ++i) {
                for (int pixel = 0; pixel < npix; ++pixel) {
                    GSkyDir dir = map.inx2dir(pixel);
                    if (regions.contains(dir)) {
                        map(pixel,i) = -1.0;
                    }
                }
            }

            // Log exclusion regions
            for (int i = 0; i < regions.size(); ++i) {
                log_value(NORMAL, "Exclusion region "+gammalib::str(i+1),
                                  region_string(*(regions[i])));
            }
        }
        else {
            log_value(NORMAL, "Exclusion regions", "None");
        }

        // Determine number of events after masking
        double sum_after = 0.0;
        for (int i = 0; i < n_ebin; ++i) {
            for (int pixel = 0; pixel < npix; ++pixel) {
                if (map(pixel, i) >= 0.0) {
                    sum_after += map(pixel, i);
                }
            }
        }

        // Write number of events before and after masking into logger
        log_value(NORMAL, "Events before masking", sum_before);
        log_value(NORMAL, "Events after masking", sum_after);

        // Put back map into the event cube
        cube->counts(map);

    } // endif: observation contained an event cube

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return region string
 *
 * @param[in] region Sky region.
 *
 * Returns a formatted region string for logging.
 ***************************************************************************/
std::string ctcubemask::region_string(const GSkyRegion& region) const
{
    // Initialise region string
    std::string rstring;

    // Get pointer to circle
    const GSkyRegionCircle* circle = dynamic_cast<const GSkyRegionCircle*>(&region);

    // If the sky region is a circle then append Right Ascenscion, Declination
    // and Radius to the region string
    if (circle != NULL) {
        rstring.append("RA="+gammalib::str(circle->centre().ra_deg())+" deg, "+
                       "DEC="+gammalib::str(circle->centre().dec_deg())+" deg, "+
                       "Radius="+gammalib::str(circle->radius())+" deg");
    }

    // Return region string
    return rstring;
}


/***********************************************************************//**
 * @brief Set output file name.
 *
 * @param[in] filename Input file name.
 *
 * This converts an input filename into an output filename by prepending a
 * prefix to the input filename. Any path will be stripped from the input
 * filename.
 ***************************************************************************/
std::string ctcubemask::set_outfile_name(const std::string& filename) const
{
    // Split input filename into path elements
    std::vector<std::string> elements = gammalib::split(filename, "/");

    // The last path element is the filename
    std::string outname = m_prefix + elements[elements.size()-1];
    
    // Return output filename
    return outname;
}


/***********************************************************************//**
 * @brief Save counts cube in FITS format.
 *
 * Save the counts cube as a FITS file. The filename of the FITS file is
 * specified by the m_outfile member.
 ***************************************************************************/
void ctcubemask::save_fits(void)
{
    // Save only if filename is non-empty and if there are observations
    if (!m_outcube.is_empty() && m_obs.size() > 0) {

        // Get CTA observation from observation container
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

        // Handle only CTA observations
        if (obs != NULL) {

            // Log counts cube file name
            log_value(NORMAL, "Counts cube file", m_outcube.url());

            // Save counts cube
            obs->save(m_outcube, clobber());
        
        } // endif: observation was a CTA observation

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save counts cube(s) in XML format.
 *
 * Save the counts cube(s) into FITS files and write the file path
 * information into an XML file. The filename of the XML file is specified by
 * the m_outfile member, the filename(s) of the counts cube(s) are built by
 * prepending the prefix given by the m_prefix member to the input counts
 * cube(s) filenames. Any path present in the input filename will be
 * stripped, i.e. the counts cube(s) will be written in the local working
 * directory (unless a path is specified in the m_prefix member).
 ***************************************************************************/
void ctcubemask::save_xml(void)
{
    // Issue warning if output filename has no .xml suffix
    log_string(TERSE, warn_xml_suffix(m_outcube.url()));

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Handle only CTA observations
        if (obs != NULL) {

            // Set event output file name
            std::string outfile = set_outfile_name(m_infiles[i]);

            // Store output file name in observation
            obs->eventfile(outfile);

            // Log counts cube file name
            log_value(NORMAL, "Counts cube file", outfile);

            // Save counts cube
            obs->save(outfile, clobber());

        } // endif: observation was a CTA observations

    } // endfor: looped over observations

    // Save observations in XML file
    m_obs.save(m_outcube);

    // Return
    return;
}
