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
#define G_GET_PARAMETERS                       "ctcubemask::get_parameters()"

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
ctcubemask::ctcubemask(void) : ctool(CTCUBEMASK_NAME, CTCUBEMASK_VERSION)
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
            ctool(CTCUBEMASK_NAME, CTCUBEMASK_VERSION)
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
ctcubemask::ctcubemask(int argc, char *argv[]) : 
            ctool(CTCUBEMASK_NAME, CTCUBEMASK_VERSION, argc, argv)
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
ctcubemask::ctcubemask(const ctcubemask& app) : ctool(app)
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
void ctcubemask::clear(void)
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

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        log.header1("Observations");
        log << m_obs << std::endl;
    }

    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Apply mask");
    }

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

            // Write header for observation
            if (logTerse()) {
                if (obs->name().length() > 1) {
                    log.header3("Observation "+obs->name());
                }
                else {
                    log.header3("Observation");
                }
            }

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

    // Write observation(s) into logger
    if (logTerse()) {
        log << std::endl;
        log.header1("Observations after event bin masking");
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
    // Write header
    if (logTerse()) {
        log << std::endl;
        if (m_obs.size() > 1) {
            log.header1("Save counts cubes");
        }
        else {
            log.header1("Save counts cube");
        }
    }

    // Get counts cube filename
    m_outcube = (*this)["outcube"].filename();

    // Case A: Save event file(s) and XML metadata information
    if (m_use_xml) {

        // Get prefix
        m_prefix = (*this)["prefix"].string();

        // Save XML
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
 * @brief Publish counts cube
 *
 * @param[in] name Counts cube name.
 *
 * Publishes the first counts cube in the observation container on the VO
 * Hub.
 ***************************************************************************/
void ctcubemask::publish(const std::string& name)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Publish counts cube");
    }

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

            // Log filename
            if (logTerse()) {
                log << "Publish \""+user_name+"\" counts cube." << std::endl;
            }

            // Publish counts cube
            cube->counts().publish(user_name);

        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 *  @exception GException::invalid_value
 *          Parameter "inobs" is required for ctcubemask.
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

    // If there are no observations in container then load them via user
    // parameters
    if (m_obs.size() == 0) {

        // Throw exception if no input observation file is given
        require_inobs(G_GET_PARAMETERS);

        // Throw exception if event list is given
        require_inobs_nolist(G_GET_PARAMETERS);

        // Build observation container without response (not needed)
        m_obs = get_observations(false);

    } // endif: there was no observation in the container

    // Get parameters
	m_regfile = (*this)["regfile"].filename();
    m_usepnt  = (*this)["usepnt"].boolean();
    if (!m_usepnt) {

        // Check RA/DEC parameters for validity to read
        if ((*this)["ra"].is_valid() && (*this)["dec"].is_valid()) {
           m_ra         = (*this)["ra"].real();
           m_dec        = (*this)["dec"].real();
           m_select_roi = true;
        }
        else {
           m_select_roi = false;
        }
    }

    // Check if radius is valid for a RoI selection
    if (m_select_roi && (*this)["rad"].is_valid()) {
       m_rad = (*this)["rad"].real();
    }
    else {
       m_select_roi = false;
    }

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
    GCTAEventCube* cube = dynamic_cast<GCTAEventCube*>(const_cast<GEvents*>(obs->events()));
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
        if (logTerse()) {
            log << gammalib::parformat("Selected energy band");
            log << ebounds.emin(e_idx1).TeV() << " - ";
            log << ebounds.emax(e_idx2).TeV() << " TeV";
            log << std::endl;
        }

        // Set all pixels inside selected energy bands but outside ROI
        // to -1.0 if requested
        if (m_select_roi) {
            GSkyRegionCircle roi(m_ra, m_dec, m_rad);
            for (int i = e_idx1; i <= e_idx2; ++i) {
                for (int pixel = 0; pixel < npix; ++pixel) {
                    GSkyDir dir = map.inx2dir(pixel);
                    if (!roi.contains(dir)) {
                        map(pixel,i) = -1.0;
                    }
                }
            }

            // Log selected energy band
            if (logTerse()) {
                log << gammalib::parformat("Selected ROI");
                log << "RA=" << m_ra << " deg, ";
                log << "DEC=" << m_dec << " deg, ";
                log << "Radius=" << m_rad << " deg";
                log << std::endl;
            }
        }

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
            if (logTerse()) {
                log << gammalib::parformat("Exclusion regions");
                log << std::endl;
                log.indent(1);
                log << regions.print();
                log << std::endl;
                log.indent(0);
            }
        }
        else {
            if (logTerse()) {
                log << gammalib::parformat("Exclusion regions");
                log << "None" << std::endl;
            }
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

        // Dump number of events before and after masking
        if (logTerse()) {
            log << gammalib::parformat("Events before masking");
            log << sum_before << std::endl;
            log << gammalib::parformat("Events after masking");
            log << sum_after << std::endl;
        }

        // Put back map into the event cube
        cube->counts(map);

    } // endif: observation contained an event cube

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
    m_ra      = 0.0;
    m_dec     = 0.0;
    m_rad     = 0.0;
    m_emin    = 0.0;
    m_emax    = 0.0;
    m_publish = false;

    // Initialise protected members
    m_obs.clear();
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
    m_ra      = app.m_ra;
    m_dec     = app.m_dec;
    m_rad     = app.m_rad;
    m_emin    = app.m_emin;
    m_emax    = app.m_emax;
    m_publish = app.m_publish;

    // Copy protected members
    m_obs           = app.m_obs;
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
    // Save only if filename is non-empty
    if (!m_outcube.is_empty()) {

        // Get CTA observation from observation container
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

        // Log filename
        if (logTerse()) {
            log << gammalib::parformat("Counts cube file");
            log << m_outcube.url() << std::endl;
        }

        // Save event list
        save_counts_map(obs, m_outcube);
        
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save counts map(s) in XML format.
 *
 * Save the counts map(s) into FITS files and write the file path information
 * into a XML file. The filename of the XML file is specified by the
 * m_outfile member, the filename(s) of the counts map(s) are built by
 * prepending the prefix given by the m_prefix member to the input counts
 * map(s) filenames. Any path present in the input filename will be stripped,
 * i.e. the counts map(s) will be written in the local working directory
 * (unless a path is specified in the m_prefix member).
 ***************************************************************************/
void ctcubemask::save_xml(void)
{
    // Issue warning if output filename has no .xml suffix
    std::string suffix = gammalib::tolower(m_outcube.url().substr(m_outcube.length()-4,4));
    if (suffix != ".xml") {
        log << "*** WARNING: Name of observation definition output file \""+
               m_outcube+"\"" << std::endl;
        log << "*** WARNING: does not terminate with \".xml\"." << std::endl;
        log << "*** WARNING: This is not an error, but might be misleading."
               " It is recommended" << std::endl;
        log << "*** WARNING: to use the suffix \".xml\" for observation"
               " definition files." << std::endl;
    }

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

            // Log filename
            if (logTerse()) {
                log << gammalib::parformat("Counts cube file");
                log << outfile << std::endl;
            }

            // Save event list
            save_counts_map(obs, outfile);

        } // endif: observation was a CTA observations

    } // endfor: looped over observations

    // Save observations in XML file
    m_obs.save(m_outcube);

    // Return
    return;
}



/***********************************************************************//**
 * @brief Save a single counts map into a FITS file
 *
 * @param[in] obs Pointer to CTA observation.
 * @param[in] outfile Output file name.
 *
 * This method saves a single counts map into a FITS file. The method does
 * nothing if the observation pointer is not valid.
 ***************************************************************************/
void ctcubemask::save_counts_map(const GCTAObservation* obs,
                                 const std::string&     outfile) const
{
    // Save only if observation is valid
    if (obs != NULL) {

        // Save observation into FITS file
        obs->save(outfile, clobber());

    } // endif: observation was valid

    // Return
    return;
}

