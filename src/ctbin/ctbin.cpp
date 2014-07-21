/***************************************************************************
 *                      ctbin - CTA data binning tool                      *
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
 * @brief CTA data binning tool implementation
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
#define G_BIN_EVENTS                    "ctbin::bin_events(GCTAObservation*)"
#define G_GET_EBOUNDS                                  "ctbin::get_ebounds()"
#define G_INIT_CUBE                                      "ctbin::init_cube()"
#define G_FILL_CUBE                      "ctbin::fill_cube(GCTAObservation*)"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_MERGE_EVENTS


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
 * param[in] obs Observation container.
 *
 * This method creates an instance of the class by copying an existing
 * observations container.
 ***************************************************************************/
ctbin::ctbin(const GObservations& obs) :
       GApplication(CTBIN_NAME, CTBIN_VERSION)
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

    // Initialise observation counter
    int n_observations = 0;

    // Compile option: if cube merging is requested, initialise the cube
    #if defined(G_MERGE_EVENTS)
    init_cube();
    #endif

    // Loop over all observations in the container
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

            // Increment number of observations
            n_observations++;

            // Save event file name (for possible saving)
            m_infiles[i] = obs->eventfile();

            // Compile option: if cube merging is requested, fill the cube;
            // otherwise just bin the events
            #if defined(G_MERGE_EVENTS)
            fill_cube(obs);
            #else
            bin_events(obs);
            #endif

        } // endif: CTA observation found

    } // endfor: looped over observations

    // Compile option: if cube merging is requested, set a single cube in
    // the observation container
    #if defined(G_MERGE_EVENTS)
    obs_cube();
    n_observations = 1;
    #endif

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

    // Compile option: if merging was requested, always save into a single
    // FITS file
    #if defined(G_MERGE_EVENTS)
    save_fits();
    #else

    /*
    // Case A: Save counts map(s) and XML metadata information
    if (m_use_xml) {
        save_xml();
    }
    // Case B: Save counts map as FITS file
    else {
        save_fits();
    }
    */
    
    // Always save FITS and an XML summary
    save_xml();
    #endif

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
            obs.load(m_evfile);

            // Append CTA observation to container
            m_obs.append(obs);

            // Signal that no XML file should be used for storage
            m_use_xml = false;

        }

        // ... otherwise try to open as XML file
        catch (GException::fits_open_error &e) {

            // Load observations from XML file
            m_obs.load(m_evfile);

            // Signal that XML file should be used for storage
            m_use_xml = true;

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
    m_ebinalg  = (*this)["ebinalg"].string();

    // If we have the binning given by a file then read filename
    if (m_ebinalg == "FILE") {
    	m_ebinfile = (*this)["ebinfile"].filename();
    }
    
    // ... otherwise read emin, emax and nebins
    else{
    	m_emin     = (*this)["emin"].real();
    	m_emax     = (*this)["emax"].real();
    	m_enumbins = (*this)["enumbins"].integer();
    }

    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (m_read_ahead) {
        m_outfile = (*this)["outfile"].filename();
        m_prefix  = (*this)["prefix"].string();
    }

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

        // Initialise empty Ebounds object
        GEbounds ebds;

        // First check wether we have FITS-file defining
        // the binning
        if (m_ebinalg == "FILE") {

        	// Open fits file to check which extension is given
        	GFits file(m_ebinfile);

        	if (file.contains("EBOUNDS")) {

        		// Close file and load ebounds
        		file.close();
        		ebds.load(m_ebinfile,"EBOUNDS");

        	} // endif: EBOUNDS extension was given

        	else if (file.contains("ENERGYBINS")) {

        		// Close file and load ebounds
        		file.close();
        		ebds.load(m_ebinfile,"ENERGYBINS");

        	} // endelse: ENERGYBINS extension was given

        	else {

                // Close file
        		file.close();
                
                // Signal that no energy boundaries have been found
                std::string msg = "No extension with name \"EBOUNDS\" or"
                                  " \"ENERGYBINS\" found in FITS file"
                                  " \""+m_ebinfile+"\".\n"
                                  "An \"EBOUNDS\" or \"ENERGYBINS\" extension"
                                  " is required if the parameter \"ebinalg\""
                                  " is set to \"FILE\".";
                throw GException::invalid_value(G_BIN_EVENTS, msg);

        	} // endelse: no suited extension found

        	// Set enumbins parameter to number of ebounds
        	m_enumbins = ebds.size();

        } // endif: ebinalg was "FILE"

        else {

        	// Initialise log mode for ebinning
        	bool log = true;

        	// check if algorithm is linear
        	if (m_ebinalg == "LIN") {
        		log = false;
        	}

        	// todo: should we also check if m_ebinalg is "LOG"
        	// and throw an exception if neither LIN/LOG/FILE
        	// is given?

        	// Setup energy range covered by data
			GEnergy  emin(m_emin, "TeV");
			GEnergy  emax(m_emax, "TeV");
        	GEbounds ebounds(m_enumbins, emin, emax,log);
        	ebds = ebounds;

        } //endif: ebinalg was not "FILE"

        // Get Good Time intervals
        GGti gti = obs->events()->gti();

        // Get map centre
        double xref = m_xref;
        double yref = m_yref;
        if (m_usepnt) {

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

        } // endif: used pointing

        // Create skymap
        GSkymap map = GSkymap(m_proj, m_coordsys,
                              xref, yref, -m_binsz, m_binsz,
                              m_nxpix, m_nypix, m_enumbins);

        // Initialise binning statistics
        int num_outside_map  = 0;
        int num_outside_ebds = 0;
        int num_in_map       = 0;

        // Fill sky map
        GCTAEventList* events =
            static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));
        for (int i = 0; i < events->size(); ++i) {

            // Get event
            GCTAEventAtom* event = (*events)[i];

            // Determine sky pixel
            GCTAInstDir* inst  = (GCTAInstDir*)&(event->dir());
            GSkyDir      dir   = inst->dir();
            GSkyPixel    pixel = map.dir2pix(dir);

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
            log << gammalib::parformat("Events in list");
            log << obs->events()->size() << std::endl;
            log << gammalib::parformat("Events in map");
            log << num_in_map << std::endl;
            log << gammalib::parformat("Events outside map area");
            log << num_outside_map << std::endl;
            log << gammalib::parformat("Events outside energy bins");
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
        obs->events(cube);

    } // endif: observation was valid

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
    m_ebounds.clear();
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
                     m_nxpix, m_nypix, m_enumbins);

    // Set energy boundaries
    get_ebounds();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill events into counts cube
 *
 * @param[in] obs CTA observation.
 *
 * @exception GException::no_list
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
        if (dynamic_cast<const GCTAEventList*>(obs->events()) == NULL) {
            throw GException::no_list(G_FILL_CUBE);
        }

        // Initialise binning statistics
        int num_outside_map  = 0;
        int num_outside_ebds = 0;
        int num_in_map       = 0;

        // Fill sky map
        GCTAEventList* events =
            static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));
        for (int i = 0; i < events->size(); ++i) {

            // Get event
            GCTAEventAtom* event = (*events)[i];

            // Determine sky pixel
            GCTAInstDir* inst  = (GCTAInstDir*)&(event->dir());
            GSkyDir      dir   = inst->dir();
            GSkyPixel    pixel = m_cube.dir2pix(dir);

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
 * @brief Put cube as single element in observation container.
 ***************************************************************************/
void ctbin::obs_cube(void)
{
    // Clear observation container
    m_obs.clear();

    // Allocate CTA observation. We need this to make sure that all
    // attributes are set correctly
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

    // Append observation to container
    m_obs.append(obs);

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
    m_ebinalg.clear();
    m_ebinfile.clear();
    m_usepnt   = false;
    m_emin     = 0.0;
    m_emax     = 0.0;
    m_enumbins = 0;
    m_xref     = 0.0;
    m_yref     = 0.0;
    m_binsz    = 0.0;
    m_nxpix    = 0;
    m_nypix    = 0;

    // Initialise protected members
    m_obs.clear();
    m_infiles.clear();
    m_use_xml    = false;
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
    m_prefix   = app.m_prefix;
    m_proj     = app.m_proj;
    m_coordsys = app.m_coordsys;
    m_ebinalg  = app.m_ebinalg;
    m_ebinfile = app.m_ebinfile;
    m_usepnt   = app.m_usepnt;
    m_emin     = app.m_emin;
    m_emax     = app.m_emax;
    m_enumbins = app.m_enumbins;
    m_xref     = app.m_xref;
    m_yref     = app.m_yref;
    m_binsz    = app.m_binsz;
    m_nxpix    = app.m_nxpix;
    m_nypix    = app.m_nypix;

    // Copy protected members
    m_obs        = app.m_obs;
    m_infiles    = app.m_infiles;
    m_use_xml    = app.m_use_xml;
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
    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get the energy boundaries
 *
 * @exception GException::invalid_value
 *            No valid energy boundary extension found.
 *
 * Get the energy boundaries according to the user parameters. The method
 * supports loading of energy boundary information from the EBOUNDS or
 * ENERGYBINS extension, or setting energy boundaries using a linear
 * or logarithmical spacing.
 *
 * @todo Ultimately this method should go in a support library as it is
 * used by several ctools.
 ***************************************************************************/
void ctbin::get_ebounds(void)
{
    // Initialse energy boundary information
    m_ebounds.clear();

    // Check whether energy binning information should be read from a FITS
    // file ...
    if (m_ebinalg == "FILE") {

        // Open energy boundary file using the EBOUNDS or ENERGYBINS
        // extension. Throw an exception if opening fails.
        GFits file(m_ebinfile);
        if (file.contains("EBOUNDS")) {
            file.close();
            m_ebounds.load(m_ebinfile,"EBOUNDS");
        }
        else if (file.contains("ENERGYBINS")) {
            file.close();
            m_ebounds.load(m_ebinfile,"ENERGYBINS");
        }
        else {
            file.close();
            std::string msg = "No extension with name \"EBOUNDS\" or"
                              " \"ENERGYBINS\" found in FITS file"
                              " \""+m_ebinfile+"\".\n"
                              "An \"EBOUNDS\" or \"ENERGYBINS\" extension"
                              " is required if the parameter \"ebinalg\""
                              " is set to \"FILE\".";
            throw GException::invalid_value(G_BIN_EVENTS, msg);
        }
        
        // Set enumbins parameter to number of ebounds
        m_enumbins = m_ebounds.size();

    } // endif: ebinalg was "FILE"

    // ... otherwise use a linear or a logarithmically-spaced energy binning
    else {

        // Initialise log mode for ebinning
        bool log = true;

        // check if algorithm is linear
        if (m_ebinalg == "LIN") {
            log = false;
        }

        // todo: should we also check if m_ebinalg is "LOG"
        // and throw an exception if neither LIN/LOG/FILE
        // is given?

        // Setup energy range covered by data
        GEnergy  emin(m_emin, "TeV");
        GEnergy  emax(m_emax, "TeV");
        m_ebounds = GEbounds(m_enumbins, emin, emax, log);

    } //endif: ebinalg was not "FILE"

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set output file name.
 *
 * @param[in] index Index of input file.
 *
 * Converts an input filename into an output filename by prepending the
 * prefix stored in the member m_prefix to the input filename. Any path will
 * be stripped from the input filename. Also a trailing ".gz" will be
 * stripped. If no input file name exist, prefix and observation 
 * number are used.
 ***************************************************************************/
std::string ctbin::set_outfile_name(const int index) const
{
    
    // Initialize output file name
    std::string outname="";
    
    // Split input filename into path elements
    std::vector<std::string> elements = gammalib::split(m_infiles[index], "/");
    
    // If there is an observation input file name
    if (elements.size() > 0) {
        outname = m_prefix + elements[elements.size()-1];
    } 
    // ...else, use observation number
    else {
    	outname = m_prefix + m_obs[index]->id() + ".fits";
    }
    
    // Strip any ".gz"
    outname = gammalib::strip_chars(outname, ".gz");

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
    m_outfile = (*this)["outfile"].filename();

    // Get CTA observation from observation container
    GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

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
    m_outfile = (*this)["outfile"].filename();
    m_prefix  = (*this)["prefix"].string();

    // Issue warning if output filename has no .xml suffix
    std::string suffix = gammalib::tolower(m_outfile.substr(m_outfile.length()-4,4));
    if (suffix != ".xml") {
        log << "*** WARNING: Name of observation definition output file \""+
               m_outfile+"\"" << std::endl;
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
            std::string outfile = set_outfile_name(i);

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
