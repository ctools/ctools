/***************************************************************************
 *                      ctbin - CTA data binning tool                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @author J. Knoedlseder
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
 * This is the main execution method of the ctbin class. It is invoked when
 * the executable is called from command line.
 *
 * The method reads the task parameters, bins the event list(s) into counts
 * map(s), and writes the results into FITS files on disk.
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

    // Initialise observation counter
    int n_observations = 0;

    // Loop over all observations in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Initialise event input and output filenames
        m_infiles.push_back("");

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

            // Increment number of observations
            n_observations++;

            // Save event file name (for possible saving)
            m_infiles[i] = obs->eventfile();

            // Bin events into counts map
            bin_events(obs);

        } // endif: CTA observation found

    } // endfor: looped over observations

    // If more than a single observation has been handled then make sure
    // that an XML file will be used for storage
    if (n_observations > 1) {
        m_use_xml = true;
    }

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
 * @brief Save counts map(s)
 *
 * This method saves the counts map(s) into FITS file(s). There are two
 * modes, depending on the m_use_xml flag.
 *
 * If m_use_xml is true, all counts map(s) will be saved into FITS files,
 * where the output filenames are constructued from the input filenames by
 * prepending the m_prefix string to name. Any path information will be
 * stripped form the input name, hence event files will be written into the
 * local working directory (unless some path information is present
 * in the prefix). In addition, an XML file will be created that gathers
 * the filename information for the counts map(s). If an XML file was present
 * on input, all metadata information will be copied from this input file.
 *
 * If m_use_xml is false, the counts map will be saved into a FITS file.
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

    // Case A: Save counts map(s) and XML metadata information
    if (m_use_xml) {
        save_xml();
    }

    // Case B: Save counts map as FITS file
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
void ctbin::get_parameters(void)
{
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
            obs.load_unbinned(m_evfile);

            // Append CTA observation to container
            m_obs.append(obs);

            // Signal that no XML file should be used for storage
            m_use_xml = false;
            
        }
        
        // ... otherwise try to open as XML file
        catch (GException::fits_open_error) {

            // Load observations from XML file
            m_obs.load(m_evfile);

            // Signal that XML file should be used for storage
            m_use_xml = true;

        }

        // Use the xref and yref parameters for binning (otherwise the
        // pointing direction(s) is/are used)
        m_xref = (*this)["xref"].real();
        m_yref = (*this)["yref"].real();

    } // endif: there was no observation in the container

    // Get remaining parameters
    m_emin     = (*this)["emin"].real();
    m_emax     = (*this)["emax"].real();
    m_enumbins = (*this)["enumbins"].integer();
    m_proj     = (*this)["proj"].string();
    m_coordsys = (*this)["coordsys"].string();
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
 * @exception GCTAException::no_pointing
 *            No valid CTA pointing found.
 *
 * This method bins the events found in a CTA events list into a counts map
 * and replaces the event list by the counts map in the observation. The
 * energy boundaries of the counts map are also stored in the observation's
 * energy boundary member.
 *
 * If the reference values for the map centre (m_xref, m_yref) are 9999.0,
 * the pointing direction of the observation is taken as the map centre.
 * Otherwise, the specified reference value is used.
 ***************************************************************************/
void ctbin::bin_events(GCTAObservation* obs)
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
        ebds.setlog(emin, emax, m_enumbins);

        // Get Good Time intervals
        GGti gti = obs->events()->gti();
        
        // Get map centre
        double xref;
        double yref;
        if (m_xref != 9999.0 && m_yref != 9999.0) {
            xref = m_xref;
            yref = m_yref;
        }
        else {
            
            // Get pointer on CTA pointing
            const GCTAPointing *pnt = obs->pointing();
            if (pnt == NULL) {
                throw GCTAException::no_pointing(G_BIN_EVENTS);
            }
            
            // Set reference point to pointing
            if (toupper(m_coordsys) == "GAL") {
                xref = pnt->dir().l_deg();
                yref = pnt->dir().b_deg();
            }
            else {
                xref = pnt->dir().ra_deg();
                yref = pnt->dir().dec_deg();
            }

        } // endelse: map centre set to pointing

        // Create skymap
        GSkymap map = GSkymap(m_proj, m_coordsys,
                              xref, yref, m_binsz, m_binsz,
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
    m_prefix.clear();
    m_proj.clear();
    m_coordsys.clear();
    m_emin     = 0.0;
    m_emax     = 0.0;
    m_enumbins = 0;
    m_xref     = 9999.0; // Flags unset (use pointing direction)
    m_yref     = 9999.0; // Flags unset (use pointing direction)
    m_binsz    = 0.0;
    m_nxpix    = 0;
    m_nypix    = 0;

    // Initialise protected members
    m_obs.clear();
    m_infiles.clear();
    m_use_xml = false;

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
    m_prefix   = app.m_prefix;
    m_proj     = app.m_proj;
    m_coordsys = app.m_coordsys;
    m_emin     = app.m_emin;
    m_emax     = app.m_emax;
    m_enumbins = app.m_enumbins;
    m_xref     = app.m_xref;
    m_yref     = app.m_yref;
    m_binsz    = app.m_binsz;
    m_nxpix    = app.m_nxpix;
    m_nypix    = app.m_nypix;

    // Copy protected members
    m_obs      = app.m_obs;
    m_infiles  = app.m_infiles;
    m_use_xml  = app.m_use_xml;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctbin::free_members(void)
{
    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set output file name.
 *
 * @param[in] filename Input file name.
 *
 * This method converts an input filename into an output filename by
 * prepending the prefix specified by m_prefix to the input filename. Any
 * path will be stripped from the input filename.
 ***************************************************************************/
std::string ctbin::set_outfile_name(const std::string& filename) const
{
    // Split input filename into path elements
    std::vector<std::string> elements = split(filename, "/");

    // The last path element is the filename
    std::string outname = m_prefix + elements[elements.size()-1];
    
    // Return output filename
    return outname;
}


/***********************************************************************//**
 * @brief Save counts map in FITS format.
 *
 * Save the counts map as a FITS file. The filename of the FITS file is
 * specified by the m_outfile member.
 ***************************************************************************/
void ctbin::save_fits(void)
{
    // Get output filename
    m_outfile = (*this)["outfile"].value();

    // Get CTA observation from observation container
    GCTAObservation* obs = dynamic_cast<GCTAObservation*>(&m_obs[0]);

    // Save event list
    save_counts_map(obs, m_outfile);

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
void ctbin::save_xml(void)
{
    // Get output filename and prefix
    m_outfile = (*this)["outfile"].value();
    m_prefix  = (*this)["prefix"].value();

    // Loop over all observation in the container
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get CTA observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(&m_obs[i]);

        // Handle only CTA observations
        if (obs != NULL) {

            // Set event output file name
            std::string outfile = set_outfile_name(m_infiles[i]);

            // Store output file name in observation
            obs->eventfile(outfile);

            // Save event list
            save_counts_map(obs, outfile);

        } // endif: observation was a CTA observations

    } // endfor: looped over observations

    // Save observations in XML file
    m_obs.save(m_outfile);

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
void ctbin::save_counts_map(const GCTAObservation* obs,
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
