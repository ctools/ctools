/***************************************************************************
 *                        ctool - ctool base class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Juergen Knoedlseder                         *
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
 * @file ctool.hpp
 * @brief ctool base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ctool.hpp"
#include "GTools.hpp"

/* __ Includes for memory usage determination ____________________________ */
#if defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>
#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>
#endif
#endif


/* __ Method name definitions ____________________________________________ */
#define G_GET_EBOUNDS                                  "ctool::get_ebounds()"
#define G_SET_FROM_CNTMAP              "ctool::set_from_cntmap(std::string&)"
#define G_GET_MEAN_POINTING        "ctool::get_mean_pointing(GObservations&)"
#define G_GET_CUBE                          "ctool::get_cube(GObservations&)"

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
ctool::ctool(void) : GApplication()
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Name constructor
 *
 * @param[in] name Application name.
 * @param[in] version Application version.
 ***************************************************************************/
ctool::ctool(const std::string& name, const std::string& version) :
       GApplication(name, version)
{
    // Initialise members
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Command line constructor
 *
 * @param[in] name Application name.
 * @param[in] version Application version.
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 ***************************************************************************/
ctool::ctool(const std::string& name, const std::string& version,
             int argc, char *argv[]) : 
       GApplication(name, version, argc, argv)
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
ctool::ctool(const ctool& app) : GApplication(app)
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
ctool::~ctool(void)
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
ctool& ctool::operator=(const ctool& app)
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
 * @brief Execute application
 *
 * This is the main execution method of a ctool. The method is invoked when
 * the executable is called from command line.
 ***************************************************************************/
void ctool::execute(void)
{
    // Signal that some parameters should be read ahead
    m_read_ahead = true;

    // Run the tool
    run();

    // Save the results
    save();

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
void ctool::init_members(void)
{
    // Initialise members
    m_read_ahead = false;
    m_use_xml    = false;

    // Set CTA time reference. G_CTA_MJDREF is the CTA reference MJD,
    // which is defined in GCTALib.hpp. This is somehow a kluge. We need
    // a better mechanism to implement the CTA reference MJD.
    m_cta_ref.set(G_CTA_MJDREF, "s", "TT", "LOCAL");

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
void ctool::copy_members(const ctool& app)
{
    // Copy members
    m_read_ahead = app.m_read_ahead;
    m_use_xml    = app.m_use_xml;
    m_cta_ref    = app.m_cta_ref;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctool::free_members(void)
{
    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Build a CTA observation with an empty event list
 *
 * Creates an empty CTA observation from user parameters. An empty event list
 * including RoI, GTI and energy boundary information is attached to the
 * observation. The method also sets the pointing direction using the ra and
 * dec parameter, the ROI based on ra, dec and rad, a single GTI based on
 * tmin and tmax, and a single energy boundary based on emin and emax. The
 * method furthermore sets the ontime, livetime and deadtime correction
 * factor.
 *
 * The following parameters are read:
 *      ra - Right Ascension of pointing and RoI centre (deg)
 *      dec - Declination of pointing and RoI centre (deg)
 *      rad - Radius of RoI (deg)
 *      deadc - Deadtime correction factor
 *      tmin - Start time
 *      tmax - Stop time
 *      emin - Minimum energy (TeV)
 *      emax - Maximum energy (TeV)
 ***************************************************************************/
GCTAObservation ctool::setup_obs(void)
{
    // Get CTA observation parameters
    double ra    = (*this)["ra"].real();
    double dec   = (*this)["dec"].real();
    double rad   = (*this)["rad"].real();
    double deadc = (*this)["deadc"].real();
    double t_min = (*this)["tmin"].real();
    double t_max = (*this)["tmax"].real();
    double e_min = (*this)["emin"].real();
    double e_max = (*this)["emax"].real();

    // Allocate CTA observation and empty event list
    GCTAObservation obs;
    GCTAEventList   list;

    // Set pointing direction
    GCTAPointing pnt;
    GSkyDir      skydir;
    skydir.radec_deg(ra, dec);
    pnt.dir(skydir);

    // Set ROI
    GCTAInstDir instdir(skydir);
    GCTARoi     roi(instdir, rad);

    // Set GTI
    GGti  gti(m_cta_ref);
    GTime tstart;
    GTime tstop;
    tstart.set(t_min, m_cta_ref);
    tstop.set(t_max, m_cta_ref);
    gti.append(tstart, tstop);

    // Set energy boundaries
    GEbounds ebounds;
    GEnergy  emin;
    GEnergy  emax;
    emin.TeV(e_min);
    emax.TeV(e_max);
    ebounds.append(emin, emax);

    // Set CTA event list attributes
    list.roi(roi);
    list.gti(gti);
    list.ebounds(ebounds);

    // Attach event list to CTA observation
    obs.events(list);

    // Set observation ontime, livetime and deadtime correction factor
    obs.pointing(pnt);
    obs.ontime(gti.ontime());
    obs.livetime(gti.ontime()*deadc);
    obs.deadc(deadc);

    // Return CTA observation
    return obs;

}


/***********************************************************************//**
 * @brief Build a skymap based on user parameters
 *
 * Creates an empty skymap from input user parameters.
 *
 * The following parameters are read:
 *      proj - Coordinate projection
 *      coordsys - Coordinate system
 *      binsz - Bin size (deg)
 *      nxpix - Number of pixels in Right Ascension or longitude
 *      nypix - Number of pixels in Declination or latitude
 ***************************************************************************/
GSkymap ctool::get_map()
{
    // Read task parameters
    double      xref     = (*this)["xref"].real();
    double      yref     = (*this)["yref"].real();
    std::string proj     = (*this)["proj"].string();
    std::string coordsys = (*this)["coordsys"].string();
    double      binsz    = (*this)["binsz"].real();
    int         nxpix    = (*this)["nxpix"].integer();
    int         nypix    = (*this)["nypix"].integer();

    // Initialise sky map
    GSkymap map = GSkymap(proj, coordsys, xref, yref, -binsz, binsz, 
                          nxpix, nypix, 1);

    // Return sky map
    return map;
}


/***********************************************************************//**
 * @brief Build a skymap based on user parameters
 *
 * @param[in] obs Observation container.
 *
 * Creates an empty skymap from input user parameters. The reference
 * coordinate is extracted from the mean pointing of the CTA observations
 * in the observation container (see get_mean_pointing()).
 *
 * The following parameters are read:
 *      proj - Coordinate projection
 *      coordsys - Coordinate system
 *      binsz - Bin size (deg)
 *      nxpix - Number of pixels in Right Ascension or longitude
 *      nypix - Number of pixels in Declination or latitude
 ***************************************************************************/
GSkymap ctool::get_map(const GObservations& obs)
{
    // Read task parameters
    std::string proj     = (*this)["proj"].string();
    std::string coordsys = (*this)["coordsys"].string();
    double      binsz    = (*this)["binsz"].real();
    int         nxpix    = (*this)["nxpix"].integer();
    int         nypix    = (*this)["nypix"].integer();

    // Initialise reference coordinates
    double xref = 0.0;
    double yref = 0.0;

    // Get pointing from observations
    GSkyDir pnt = get_mean_pointing(obs);

    // Use correct coordinate system
    if (gammalib::toupper(coordsys) == "GAL") {
        xref = pnt.l_deg();
        yref = pnt.b_deg();
    }
    else {
        xref = pnt.ra_deg();
        yref = pnt.dec_deg();
    }

    // Initialise sky map
    GSkymap map = GSkymap(proj, coordsys, xref, yref, -binsz, binsz, 
                          nxpix, nypix, 1);

    // Return sky map
    return map;
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
 * The following parameters are read:
 *      ebinalg - Energy binning algorithm
 *      emin - Minimum energy (if ebinalg != FILE)
 *      emax - Maximum energy (if ebinalg != FILE)
 *      enumbins - Number of energy bins (if ebinalg != FILE)
 ***************************************************************************/
GEbounds ctool::get_ebounds(void)
{
    // Allocate energy boundaries
    GEbounds ebounds;

    // Get energy binning algorithm
    std::string ebinalg = (*this)["ebinalg"].string();

    // If energy binning algorithm is of type "FILE" (case sensitive), then
    // read energy boundaries from FITS file ...
    if (ebinalg == "FILE") {

        // Get filename
        std::string ebinfile = (*this)["ebinfile"].filename();

        // Open energy boundary file using the EBOUNDS or ENERGYBINS
        // extension. Throw an exception if opening fails.
        GFits file(ebinfile);
        if (file.contains("EBOUNDS")) {
            file.close();
            ebounds.load(ebinfile,"EBOUNDS");
        }
        else if (file.contains("ENERGYBINS")) {
            file.close();
            ebounds.load(ebinfile,"ENERGYBINS");
        }
        else {
            file.close();
            std::string msg = "No extension with name \"EBOUNDS\" or"
                              " \"ENERGYBINS\" found in FITS file"
                              " \""+ebinfile+"\".\n"
                              "An \"EBOUNDS\" or \"ENERGYBINS\" extension"
                              " is required if the parameter \"ebinalg\""
                              " is set to \"FILE\".";
            throw GException::invalid_value(G_GET_EBOUNDS, msg);
        }
        
    } // endif: ebinalg was "FILE"

    // ... otherwise use a linear or a logarithmically-spaced energy binning
    else {

        // Get task parameters
    	double emin     = (*this)["emin"].real();
    	double emax     = (*this)["emax"].real();
    	int    enumbins = (*this)["enumbins"].integer();

        // Initialise log mode for ebinning
        bool log = true;

        // Check if algorithm is linear
        if (ebinalg == "LIN") {
            log = false;
        }

        // Setup energy bins
        ebounds = GEbounds(enumbins, GEnergy(emin, "TeV"),
                                     GEnergy(emax, "TeV"), log);

    } // endelse: ebinalg was not "FILE"

    // Return energy boundaries
    return ebounds;
}


/***********************************************************************//**
 * @brief Get observation container
 *
 * @param[in] get_response Indicates whether response information should
 *                         be loaded (default: true)
 *
 * Get an observation container according to the user parameters. The method
 * supports loading of a individual FITS file or an observation definition
 * file in XML format. If the input filename is empty, parameters are read
 * to build a CTA observation from scratch.
 ***************************************************************************/
GObservations ctool::get_observations(const bool& get_response)
{
    // Initialise empty observation container
    GObservations obs;

    // Get the filename from the input parameters
    std::string filename = (*this)["inobs"].filename();

    // If no observation definition file has been specified then read all
    // parameters that are necessary to create an observation from scratch
    // (see method setup_obs)
    if ((filename == "NONE") || (gammalib::strip_whitespace(filename) == "")) {

        // Setup a new CTA observation
        GCTAObservation cta_obs = setup_obs();

        // Get response if required
        if (get_response) {

            // Set response
           set_obs_response(&cta_obs);

        } // endif: response was required

       // Append observation to container
       obs.append(cta_obs);

    } // endif: filename was "NONE" or ""

    // ... otherwise we have a file name
    else {

        // Try first to open as FITS file
        try {

            // Allocate empyt CTA observation
            GCTAObservation cta_obs;

            // Load data
            cta_obs.load(filename);

            // Get response if required
            if (get_response) {

                // Set response
                set_obs_response(&cta_obs);

            } // endif: response was required

            // Append observation to container
            obs.append(cta_obs);

            // Signal that no XML file should be used for storage
            m_use_xml = false;

        }

        // ... otherwise try to open as XML file
        catch (GException::fits_open_error &e) {

            // Load observations from XML file
            obs.load(filename);

            // Get response if required
            if (get_response) {

                // For all observations that have no response, set the response
                // from the task parameters
                set_response(obs);

            } // endif: response was required

            // Signal that XML file should be used for storage
            m_use_xml = true;

        } // endcatch: file was an XML file

    }

    // Return observation container
    return obs;

}


/***********************************************************************//**
 * @brief Derives mean pointing from CTA observations
 *
 * @param[in] obs Observation container.
 * @return Pointing direction.
 *
 * Computes the mean pointing direction from all CTA observations in the
 * observation container.
 *
 * @todo This method does not work properly for wrap arounds in Right
 *       Ascension. It certainly will also lead to bad results near the
 *       celestial pole.
 ***************************************************************************/
GSkyDir ctool::get_mean_pointing(const GObservations& obs)
{
    // Initialise values
    double xref  = 0.0;
    double yref  = 0.0;
    int    n_pnt = 0;

    // Loop over all observations to compute mean pointings
    for (int i = 0; i < obs.size(); ++i) {

        // If we have a CTA observation then extract the pointing direction
        const GCTAObservation* cta_obs = dynamic_cast<const GCTAObservation*>(obs[i]);
        if (cta_obs != NULL) {

           // Get pointing
           const GCTAPointing& pnt = cta_obs->pointing();

           // Add to coordinates
           xref  += pnt.dir().ra_deg();
           yref  += pnt.dir().dec_deg();
           n_pnt += 1;

       } // endif: cta observation was valid

    } // endfor: loop over all observations

    // Signal if no pointing is found
    if (n_pnt < 1) {
        std::string msg = "No valid CTA observation has been found in "
                          "observation list, hence no pointing information "
                          "could be extracted. Use the \"usepnt=no\" "
                          "option and specify the pointing explicitly.";
        throw GException::invalid_value(G_GET_MEAN_POINTING, msg);
    }

    // Calculate mean coordinates
    double ra  = xref / double(n_pnt);
    double dec = yref / double(n_pnt);

    // Initialise sky dir
    GSkyDir dir = GSkyDir();
    dir.radec_deg(ra,dec);

    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief Build CTA event cube based on user parameters
 *
 * Creates an empty CTA event cube from input user parameters. It appends a
 * dummy GTI to the event cube.
 *
 * The following parameters are read:
 *      xref - Reference in Right Ascension or longitude (deg)
 *      yref - Reference in Declination or latitude (deg)
 *      proj - Coordinate projection
 *      coordsys - Coordinate system
 *      binsz - Bin size (deg)
 *      nxpix - Number of pixels in Right Ascension or longitude
 *      nypix - Number of pixels in Declination or latitude
 *      ebinalg - Energy binning algorithm
 *      emin - Minimum energy (if ebinalg != FILE)
 *      emax - Maximum energy (if ebinalg != FILE)
 *      enumbins - Number of energy bins (if ebinalg != FILE)
 ***************************************************************************/
GCTAEventCube ctool::get_cube(void)
{
    // Get skymap
    GSkymap map = get_map();

    // Set energy boundaries
    GEbounds ebounds  = get_ebounds();

    // Set dummy GTI that is needed for event cube creation
    GGti gti;
    gti.append(GTime(0.0), GTime(0.1234));

    // Extend skymap to the requested number of maps
    map.nmaps(ebounds.size());

    // Build event cube from given information
    GCTAEventCube cube = GCTAEventCube(map, ebounds, gti);

    // Return cube
    return cube;
}


/***********************************************************************//**
 * @brief Build CTA event cube based on user parameters
 *
 * @param[in] obs Observation container.
 *
 * Creates an empty CTA event cube from input user parameters. The reference
 * coordinate is extracted from the mean pointing of the CTA observations
 * in the observation container (see get_mean_pointing()). This method will
 * be used if a ctool has the input parameter usepnt=true. The method
 * appends a dummy GTI to the event cube.
 *
 * The following parameters are read:
 *      proj - Coordinate projection
 *      coordsys - Coordinate system
 *      binsz - Bin size (deg)
 *      nxpix - Number of pixels in Right Ascension or longitude
 *      nypix - Number of pixels in Declination or latitude
 *      ebinalg - Energy binning algorithm
 *      emin - Minimum energy (if ebinalg != FILE)
 *      emax - Maximum energy (if ebinalg != FILE)
 *      enumbins - Number of energy bins (if ebinalg != FILE)
 ***************************************************************************/
GCTAEventCube ctool::get_cube(const GObservations& obs)
{
    // Get skymap
    GSkymap map = get_map(obs);

    // Set energy boundaries
    GEbounds ebounds = get_ebounds();

    // Set dummy GTI that is needed for event cube creation
    GGti gti;
    gti.append(GTime(0.0), GTime(0.1234));

    // Extend skymap to the requested number of maps
    map.nmaps(ebounds.size());

    // Initialise event cube from all parameters
    GCTAEventCube cube = GCTAEventCube(map, ebounds, gti);

    // Return event cube
    return cube;
}


/***********************************************************************//**
 * @brief Get current resident set size (physical memory use) in Bytes
 *
 * @return Physical memory use in Bytes.
 ***************************************************************************/
size_t ctool::get_current_rss(void)
{
    // Initialize resident set size
    size_t rss = 0;

    // Determine resident set size (architecture dependent)
    // OSX
    #if defined(__APPLE__) && defined(__MACH__)
    #ifdef MACH_TASK_BASIC_INFO
    struct mach_task_basic_info info;
    mach_msg_type_number_t      infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount) == KERN_SUCCESS) {
        rss = (size_t)info.resident_size;
    }
    #else
    struct task_basic_info info;
    mach_msg_type_number_t info_count = TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), TASK_BASIC_INFO,
        (task_info_t)&info, &info_count) == KERN_SUCCESS) {
        rss = (size_t)info.resident_size;
    }
    #endif
    // Linux
    #elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    FILE* fp = NULL;
    if ((fp = fopen( "/proc/self/statm", "r" )) != NULL) {
        if (fscanf( fp, "%*s%ld", &rss ) == 1) {
            rss *= (size_t)sysconf(_SC_PAGESIZE);
        }
        fclose(fp);
    }
    // AIX, BSD, Solaris, and Unknown OS
    #else
    rss = 0;
    #endif

    // Return resident set size
    return rss;
}


/***********************************************************************//**
 * @brief Get a CTA ebent cube from from file
 *
 * @param[in] filename FITS file of count cube.
 *
 * @exception GException::invalid_value
 *            Counts map not of WCS type
 *
 * Gets an event cube from a FITS file. This function is used in many ctools
 * to obtain a skymap definition from an event cube.
 *
 * @todo Is it necessary to restrict this method to WCS coordinates?
 ***************************************************************************/
GCTAEventCube ctool::set_from_cntmap(const std::string filename)
{
    // Allocate CTA observation
    GCTAEventCube cube;

    // Load counts map in CTA observation
    cube.load(filename);

    // Check for valid WCS coordinate system
    const GWcs* wcs = dynamic_cast<const GWcs*>(cube.map().projection());
    if (wcs == NULL) {
        std::string msg = "Counts map project is not of WCS type.";
        throw GException::invalid_value(G_SET_FROM_CNTMAP, msg);
    }

    // Return cube
    return cube;
}


/***********************************************************************//**
 * @brief Set response for all CTA observations in container
 *
 * @param[in,out] obs Observation container
 *
 * Set the response for a CTA observations in the container that so far have
 * no response using the "caldb" and "irf" task parameters.
 *
 * The following parameters are read:
 *      caldb - Calibration database
 *      irf - Instrument response function
 ***************************************************************************/
void ctool::set_response(GObservations& obs)
{
    // Setup response for all observations
    for (int i = 0; i < obs.size(); ++i) {

        // Is this observation a CTA observation?
        GCTAObservation* cta = dynamic_cast<GCTAObservation*>(obs[i]);

        // Yes ...
        if (cta != NULL) {

            // Set response if we don't have one
            if (!cta->has_response()) {
                set_obs_response(cta);
            }

        } // endif: observation was a CTA observation

    } // endfor: looped over all observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set response for CTA observation
 *
 * @param[in,out] obs CTA observation
 *
 * Set the response for a CTA observation using the "caldb" and "irf"
 * task parameters.
 *
 * The following parameters are read:
 *      caldb - Calibration database
 *      irf - Instrument response function
 ***************************************************************************/
void ctool::set_obs_response(GCTAObservation* obs)
{
    // Load response information
    std::string database = (*this)["caldb"].string();
    std::string irf      = (*this)["irf"].string();

    // Set calibration database. If "database" is a valid directory then use
    // this as the pathname to the calibration database. Otherwise, interpret
    // "database" as the instrument name, the mission being "cta"
    GCaldb caldb;
    if (gammalib::dir_exists(database)) {
        caldb.rootdir(database);
    }
    else {
        caldb.open("cta", database);
    }

    // Set reponse
    obs->response(irf, caldb);

    // Return
    return;
}
