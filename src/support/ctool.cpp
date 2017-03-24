/***************************************************************************
 *                        ctool - ctool base class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2017 by Juergen Knoedlseder                         *
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
#include <cstdlib>         // std::getenv() function
#include <cstdio>          // std::fopen(), etc. functions
#include <clocale>         // std::setlocale function
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
#define G_SETUP_OBSERVATION       "ctool::setup_observations(GObservations&)"
#define G_SETUP_MODELS    "ctool::setup_models(GObservations&, std::string&)"
#define G_GET_MEAN_POINTING        "ctool::get_mean_pointing(GObservations&)"
#define G_CREATE_EBOUNDS                            "ctool::create_ebounds()"
#define G_RESTORE_EDISP               "ctool::restore_edisp(GObservations&, "\
                                                        "std::vector<bool>&)"
#define G_SET_OBS_RESPONSE        "ctool::set_obs_response(GCTAObservation*)"
#define G_PROVIDE_HELP                                "ctool::provide_help()"

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Name constructor
 *
 * @param[in] name Application name.
 * @param[in] version Application version.
 *
 * Constructs a ctool from the @p name and @p version. The constructor uses
 * the equivalent GApplication constructor to set the parameter filename to
 * "<name>.par" and the log filename to "<name>".log. The parameters will be
 * loaded from the parameter file.
 *
 * No log file will be opened. To open the log file an explicit call to the
 * logFileOpen() method is required.
 *
 * This constructor should be used for using a ctool from Python.
 ***************************************************************************/
ctool::ctool(const std::string& name, const std::string& version) :
       GApplication(name, version)
{
    // Initialise members
    init_members();

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
 *
 * Constructs a ctool from the @p name, @p version and command line
 * arguments. The constructor uses the equivalent GApplication constructor
 * to set the parameter filename to "<name>.par" and the log filename to
 * "<name>".log. The parameters will be loaded from the parameter file.
 * In addition, the constructor opens the log file.
 *
 * If the "--help" option is provided as command line argument a help text
 * about the usage of the ctool will be shown in the console and the ctool
 * will exit. No log file will be opened in that case.
 ***************************************************************************/
ctool::ctool(const std::string& name, const std::string& version,
             int argc, char *argv[]) : 
       GApplication(name, version, argc, argv)
{
    // Catch --help option before doing anything else
    if (need_help()) {
        provide_help();
        exit(0);
    }
    
    // Initialise members
    init_members();

    // Open the log file
    logFileOpen();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Application.
 *
 * Construct a ctools from another ctool instance.
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
 *
 * Assigns one ctool to another.
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
 * the executable is called from the command line. It signals that output
 * parameters should be read ahead, runs the tools, and saves the results.
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

    // Set language to english to make sure that a dot is a dot
    std::setlocale(LC_ALL, "en_US.UTF-8");

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
    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup observation container
 *
 * @param[in] obs Observation container
 * @param[in] response Require response
 * @param[in] list Accept event list as "inobs" parameter
 * @param[in] cube Accept counts cube as "inobs" parameter
 *
 * @exception GException::invalid_value
 *            Invalid "inobs" parameter encountered.
 *
 * Setup an observation container by extracting all required missing
 * information from user parameters.
 *
 * If the observation container @p obs is empty, the method will extract the
 * filename from the "inobs" parameter and load the observation container
 * from the specified file. If the filename is either empty or "NONE" the
 * method will throw an exception.
 *
 * If the "inobs" file is a FITS file, and depending on the information that
 * is found in the FITS file, the method will append a single CTA observation
 * to the observation container, containing either an event list or a counts
 * cube. It will set the m_use_xml member to false to signal that no
 * observation XML file should be used for saving. The @p list and @p cube
 * parameters can be used to specify whether the method actually should accept
 * event lists or counts cube. If it shouldn't the method will throw an
 * exception.
 *
 * If the "inobs" file is not a FITS file, the method assumes that an
 * observation definition XML file has been specified, and will load that
 * file into the observation container. It will set the m_use_xml member to
 * true to signal that an observation XML file should be used for saving. No
 * checking on event list or counts cubes will be done on the observation
 * container.
 *
 * By default, the method will also setup the response for all CTA
 * observations in the container. In case that no response information is
 * required, the @p response argument should be set to @p false. See the
 * set_response() method for details.
 *
 *
 * The method will also setup the observation boundaries for all CTA
 * observations in the container that contain event lists. For details, see
 * the set_obs_bounds() method.
 ***************************************************************************/
void ctool::setup_observations(GObservations& obs,
                               const bool&    response,
                               const bool&    list,
                               const bool&    cube)
{
    // Load observations if there are none in the container
    if (obs.size() == 0) {

        // Get the filename from the "inobs" parameter
        GFilename filename = (*this)["inobs"].filename();

        // Throw an exception if the "inobs" parameter is not a valid
        // filename
        if (!is_valid_filename(filename)) {
            std::string msg = "The \"inobs\" parameter \""+filename.url()+
                              "\" is not a valid filename. Please specify "
                              "a valid event list, counts cube or "
                              "observation definition XML file.";
            throw GException::invalid_value(G_SETUP_OBSERVATION, msg);
        }

        // If file does not exist then throw an exceptions
        if (!filename.exists()) {
            std::string msg = "The \"inobs\" parameter \""+filename.url()+
                              "\" specifies a file that does not "
                              "exist. Please specify  a valid event "
                              "list, counts cube or  observation "
                              "definition XML file.";
            throw GException::invalid_value(G_SETUP_OBSERVATION, msg);
        }

        // If the file is a FITS file then load it into a CTA
        // observation and append that observation to the container ...
        else if (filename.is_fits()) {

            // Load CTA observation
            GCTAObservation cta(filename);

            // If no event list is accepted but an event list has been
            // loaded then throw an exception
            if (!list && cta.eventtype() == "EventList") {
                std::string msg = "The \"inobs\" parameter \""+filename.url()+
                                  "\" specifies an event list, but this "
                                  "tool does not accept event lists.";
                throw GException::invalid_value(G_SETUP_OBSERVATION, msg);
            }

            // If no counts cube is accepted but a conts cube has been
            // loaded then throw an exception
            if (!cube && cta.eventtype() == "CountsCube") {
                std::string msg = "The \"inobs\" parameter \""+filename.url()+
                                  "\" specifies a counts cube, but this "
                                  "tool does not accept counts cubes.";
                throw GException::invalid_value(G_SETUP_OBSERVATION, msg);
            }

            // Append CTA observation to observation container
            obs.append(cta);

            // Signal that no XML file should be used for storage
            m_use_xml = false;

        } // endelse: file was a FITS file

        // ... otherwise it is assumed that the file is an observation
        // definition XML file and the file is loaded into the
        // observation container
        else {

            // Load observation definition file
            obs.load(filename);

            // Signal that XML file should be used for storage
            m_use_xml = true;

        } // endelse: file was an XML file

    } // endif: there were no observations in the container

    // Setup the response for all CTA observations in the observation
    // container if response information is required
    if (response) {
        set_response(obs);
    }

    // Setup boundaries for all CTA observations containing event lists
    set_obs_bounds(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup model container
 *
 * @param[in] obs Observation container
 * @param[in] name Mandatory model name
 *
 * @exception GException::invalid_value
 *            Model definiton XML filename not valid or mandatory model name
 *            not found in model container.
 *
 * Setup a model container by loading the models from the "inmodel"
 * parameter.
 *
 * If the model container in the observation container @p obs is empty, the
 * method will extract the model container filename from the "inmodel"
 * parameter, load the model container and assign it to the observations
 * container. If the filename is either empty or "NONE" the method will
 * throw an exception.
 *
 * If a mandatory model name is specified, the method will further check
 * whether a model with the name exists in the model container. If it does
 * not exist, an exception is thrown.
 ***************************************************************************/
void ctool::setup_models(GObservations&     obs,
                         const std::string& name)
{
    // If there are no models in the observation container then load model
    // container from the "inmodel" parameter and attach it to the
    // observation container
    if (obs.models().size() == 0) {

        // Get models XML filename
        GFilename filename = (*this)["inmodel"].filename();

        // Throw an exception if the "inmodel" parameter is not a valid
        // filename
        if (!is_valid_filename(filename)) {
            std::string msg = "The \"inmodel\" parameter \""+filename.url()+
                              "\" is not a valid filename. Please specify "
                              "a valid model definition XML file.";
            throw GException::invalid_value(G_SETUP_MODELS, msg);
        }

        // Load model container and assign it to observations container
        obs.models(GModels(filename.url()));

    } // endif: no models were in observations container

    // If a mandatory model name exists then check whether this model is
    // part of the observation container. Throw an exception if the model
    // does not exist
    if (!name.empty()) {

        // Throw an exception if the model name does not exist
        if (!obs.models().contains(name)) {
            std::string msg = "Model \""+name+"\" not found in model "
                              "container. Please add a model with that name "
                              "or check for possible typos.";
            throw GException::invalid_value(G_SETUP_MODELS, msg);
        }

    } // endif: mandatory model name existed

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create energy boundaries from user parameters
 *
 * @exception GException::invalid_value
 *            No valid energy boundary extension found or invalid extension
 *            name.
 *
 * Get the energy boundaries according to the user parameters. The method
 * supports loading of energy boundary information from the `EBOUNDS` or
 * `ENERGYBINS` extension, loading of energy value information from the
 * `ENERGIES` extension, or setting energy boundaries using a linear or
 * logarithmical spacing.
 *
 * The following parameters are read:
 *
 *      ebinalg - Energy binning algorithm
 *      ebinfile - Name of file with energy binning
 *      emin - Minimum energy (if ebinalg != FILE)
 *      emax - Maximum energy (if ebinalg != FILE)
 *      enumbins - Number of energy bins (if ebinalg != FILE)
 ***************************************************************************/
GEbounds ctool::create_ebounds(void)
{
    // Allocate energy boundaries
    GEbounds ebounds;

    // Get energy binning algorithm
    std::string ebinalg = (*this)["ebinalg"].string();

    // If energy binning algorithm is of type "FILE" (case sensitive), then
    // read energy boundaries from FITS file ...
    if (ebinalg == "FILE") {

        // Get filename
        GFilename ebinfile = (*this)["ebinfile"].filename();

        // If no extension name was provided then use default extension names
        if (!ebinfile.has_extname()) {

            // Open energy boundary file using the EBOUNDS or ENERGYBINS
            // extension or energy values using the ENERGIES extension.
            // Throw an exception if opening fails.
            GFits file(ebinfile.url());
            if (file.contains(gammalib::extname_ebounds)) {
                file.close();
                ebounds.load(ebinfile);
            }
            else if (file.contains("ENERGYBINS")) {
                file.close();
                ebinfile = ebinfile.url() + "[ENERGYBINS]";
                ebounds.load(ebinfile);
            }
            else if (file.contains(gammalib::extname_energies)) {
                file.close();
                ebounds.set(GEnergies(ebinfile));
            }
            else {
                file.close();
                std::string msg = "No extension with name \""+
                                  gammalib::extname_ebounds+"\", "
                                  "\"ENERGYBINS\" or \""+
                                  gammalib::extname_energies+"\" found in FITS "
                                  "file \""+ebinfile+"\". Please specify a "
                                  "valid energy binning file.";
                throw GException::invalid_value(G_CREATE_EBOUNDS, msg);
            }
        }

        // ... otherwise load energy boundaries from filename including
        // extension
        else {

            // Open energy boundary file
            GFits file(ebinfile.url());

            // If FITS file does not contain requested extension then throw
            // an exception
            if (!file.contains(ebinfile.extname())) {
                std::string msg = "No extension \""+ebinfile.extname()+"\" "
                                  "found in energy binning file \""+
                                  ebinfile.url()+"\". Please provide a valid "
                                  "extension name.";
                throw GException::invalid_value(G_CREATE_EBOUNDS, msg);
            }

            // Get table extension
            GFitsTable& table = *file.table(ebinfile.extname());

            // If table contains one column then load energies, otherwise load
            // energy boundaries.
            if (table.ncols() == 1) {
                file.close();
                ebounds.set(GEnergies(ebinfile));
            }
            else {
                file.close();
                ebounds.load(ebinfile);
            }

        } // endelse: loaded energy boundaries from table extension
        
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
 * @brief Create a skymap from user parameters
 *
 * @param[in] obs Observation container.
 *
 * Creates an empty skymap from input user parameters. The reference
 * coordinate is extracted from the mean pointing of the CTA observations
 * in the observation container (see get_mean_pointing()).
 *
 * The following parameters are read:
 *      usepnt - Use pointing for reference coordinate?
 *      xref - Reference in Right Ascension or longitude (if usepnt=no)
 *      yref - Reference in Declination or latitude (if usepnt=no)
 *      proj - Coordinate projection
 *      coordsys - Coordinate system
 *      binsz - Bin size (deg)
 *      nxpix - Number of pixels in Right Ascension or longitude
 *      nypix - Number of pixels in Declination or latitude
 ***************************************************************************/
GSkyMap ctool::create_map(const GObservations& obs)
{
    // Read task parameters
    double xref   = 0.0;
    double yref   = 0.0;
    bool   usepnt = (*this)["usepnt"].boolean();
    if (!usepnt) {
        xref = (*this)["xref"].real();
        yref = (*this)["yref"].real();
    }
    std::string proj     = (*this)["proj"].string();
    std::string coordsys = (*this)["coordsys"].string();
    double      binsz    = (*this)["binsz"].real();
    int         nxpix    = (*this)["nxpix"].integer();
    int         nypix    = (*this)["nypix"].integer();

    // If requested, get pointing from observations
    if (usepnt) {
    
        // Get mean pointing
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

    } // endif: got mean pointing as reference

    // Initialise sky map
    GSkyMap map = GSkyMap(proj, coordsys, xref, yref, -binsz, binsz, 
                          nxpix, nypix, 1);

    // Return sky map
    return map;
}


/***********************************************************************//**
 * @brief Create a CTA event cube from user parameters
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
GCTAEventCube ctool::create_cube(const GObservations& obs)
{
    // Get skymap
    GSkyMap map = create_map(obs);

    // Set energy boundaries
    GEbounds ebounds = create_ebounds();

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
 * @brief Create a CTA observation from user parameters
 *
 * Creates an empty CTA observation from user parameters. An empty event list
 * including RoI, GTI and energy boundary information is attached to the
 * observation. The method also sets the pointing direction using the ra and
 * dec parameter, the ROI based on ra, dec and rad, a single GTI based on
 * tmin and tmax, and a single energy boundary based on emin and emax. The
 * method furthermore sets the ontime, livetime and deadtime correction
 * factor.
 *
 * The following parameters are read
 *
 *      ra:    Right Ascension of pointing and RoI centre (deg)
 *      dec:   Declination of pointing and RoI centre (deg)
 *      rad:   Radius of RoI (deg)
 *      deadc: Deadtime correction factor
 *      tmin:  Start time
 *      tmax:  Stop time
 *      emin:  Minimum energy (TeV)
 *      emax:  Maximum energy (TeV)
 ***************************************************************************/
GCTAObservation ctool::create_cta_obs(void)
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

    // Attach empty event list to CTA observation
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
 * @brief Throws exception if inobs parameter is not valid
 *
 * @param[in] method Method name.
 *
 * Throw an exception if the "inobs" parameter is either "NONE" or an empty
 * string.
 *
 * @todo This method is only used by csobsinfo, csresmap and csspec. We
 * should try to get rid of it
 ***************************************************************************/
void ctool::require_inobs(const std::string& method)
{
    // Get inobs filename
    GFilename filename = (*this)["inobs"].filename();

    // Throw exception if no infile is given
    if (!is_valid_filename(filename)) {
        std::string msg = "A valid file needs to be specified for the "
                          "\"inobs\" parameter, yet \""+filename.url()+
                          "\" was given."
                          " Specify a valid observation definition or "
                          "FITS file to proceed.";
        throw GException::invalid_value(method, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Throws exception if inobs parameter is a counts cube
 *
 * @param[in] method Method name.
 *
 * Throw an exception if the inobs parameter is a counts cube.
 *
 * @todo This method is only used by ctobssim. We should try to get rid of it
 ***************************************************************************/
void ctool::require_inobs_nocube(const std::string& method)
{
    // Get inobs filename
    GFilename filename = (*this)["inobs"].filename();

    // Continue only if we have a FITS file
    if (filename.is_fits()) {

        // Signal no cube
        bool is_cube = false;
    
        // Try loading file as counts cube. If this is successful then
        // throw an exception
        try {

            // Load cube from file
            GCTAEventCube cube(filename);

            // If we're still alive then signal that we have a cube
            is_cube = true;

        }

        // Catch any exceptions
        catch (...) {
            ;
        }

        // If we have a counts cube then throw an exception
        if (is_cube) {
            std::string msg = "A counts cube has been specified for the "
                              "\"inobs\" parameter, yet no counts cube "
                              "can be specified as input observation."
                              " Instead, specify an event list or an "
                              "observation definition file.";
            throw GException::invalid_value(method, msg);
        }
    
    } // endif: we had a FITS file
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Log observation container
 *
 * @param[in] chatter Minimum required chattiness
 * @param[in] obs Observation container
 * @param[in] what String specifying the container content
 *
 * Log observations if chattiness is at least @p chatter.
 ***************************************************************************/
void ctool::log_observations(const GChatter&      chatter,
                             const GObservations& obs,
                             const std::string&   what)
{
    // Get chattiness of ctool
    GChatter chattiness = static_cast<GChatter>((*this)["chatter"].integer());

    // Only write message if chattiness is at least equal to the minimum
    // required chattiness
    if (chattiness >= chatter) {

        // Log observation header
        log << std::endl;
        log.header1(gammalib::number(what, obs.size()));

        // Log observation content dependent on chattiness
        log << obs.print(chattiness) << std::endl;

    } // endif: Chattiness satisfied minimum required level

    // Return
    return;
}


/***********************************************************************//**
 * @brief Log model container
 *
 * @param[in] chatter Minimum required chattiness
 * @param[in] models Model container
 * @param[in] what String specifying the container content
 *
 * Log models if chattiness is at least @p chatter.
 ***************************************************************************/
void ctool::log_models(const GChatter&    chatter,
                       const GModels&     models,
                       const std::string& what)
{
    // Get chattiness of ctool
    GChatter chattiness = static_cast<GChatter>((*this)["chatter"].integer());

    // Only write message if chattiness is at least equal to the minimum
    // required chattiness
    if (chattiness >= chatter) {

        // Log observation header
        log << std::endl;
        log.header1(gammalib::number(what, models.size()));

        // Log observation content dependent on chattiness
        log << models.print(chattiness) << std::endl;

    } // endif: Chattiness satisfied minimum required level

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return RoI from user parameters
 *
 * @param[in] obs Observation container.
 * @return Region of interest.
 *
 * Returns region of interest from the user parameters "ra", "dec" and "rad".
 * If the "usepnt" parameter exists and is "yes" then postpone then assume
 * that the Right Ascension and Declination of the RoI centre will be set
 * from the pointing direction somewhere else.
 ***************************************************************************/
GCTARoi ctool::get_roi(void)
{
    // Initialise an empty RoI
    GCTARoi roi;

    // Signal whether the pointing instead of the user parameters should be
    // use for the RoI centre
    bool usepnt = (has_par("usepnt") && (*this)["usepnt"].boolean());

    // If pointing should not bet used the extract the RoI centre from the
    // "ra" and "dec" parameters
    if (!usepnt) {

        // Read "ra" and "dec" parameters only if both are valid
        if ((*this)["ra"].is_valid() && (*this)["dec"].is_valid()) {

            // Read parameters
            double ra  = (*this)["ra"].real();
            double dec = (*this)["dec"].real();

            // Set RoI
            GCTAInstDir instdir;
            instdir.dir().radec_deg(ra, dec);
            roi.centre(instdir);

            // Signal that we have a pointing direction
            usepnt = true;

        } // endif: "ra" and "dec" were valid

    } // endif: pointing was not used for RoI centre

    // If we have a RoI centre and if the "rad" parameter is valid then
    // extract the RoI radius from the parameter
    if (usepnt && (*this)["rad"].is_valid()) {
        roi.radius((*this)["rad"].real());
    }

    // Return RoI
    return roi;
}


/***********************************************************************//**
 * @brief Set response for all CTA observations in container
 *
 * @param[in,out] obs Observation container
 *
 * Set the response for all CTA observations in the container that so far
 * have no response information available. See the set_obs_response()
 * method for more information about the user parameters that are read.
 ***************************************************************************/
void ctool::set_response(GObservations& obs)
{
    // Setup response for all observations
    for (int i = 0; i < obs.size(); ++i) {

        // Is this observation a CTA observation?
        GCTAObservation* cta = dynamic_cast<GCTAObservation*>(obs[i]);

        // Yes, the set the response if the observation does not yet have one
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
 * @brief Set energy dispersion to CTA observations
 *
 * @param[in] obs Observation container.
 * @param[in] edisp Requested energy dispersion flag.
 * @return Vector of old energy dispersion flags.
 *
 * Applies energy dispersion to all CTA observations. The method returns a
 * vector with the old energy dispersion flags.
 ***************************************************************************/
std::vector<bool> ctool::set_edisp(GObservations& obs, const bool& edisp) const
{
    // Allocate vector for energy dispersion flags of all observation in
    // the container
    std::vector<bool> old_edisp(obs.size(), false);

    // Loop over all observations in observation container
    for (int i = 0; i < obs.size(); ++i) {

        // Is this observation a CTA observation?
        GCTAObservation* cta = dynamic_cast<GCTAObservation*>(obs[i]);

        // Yes, then set the energy dispersion flag
        if (cta != NULL) {
        
            // Save old energy dispersion flag
            old_edisp[i] = cta->response()->apply_edisp();

            // Set energy dispersion flag according to input parameter
            cta->response()->apply_edisp(edisp);

        } // endif: observation was CTA observation

    } // endfor: loop over all observations in container

    // Return old energy dispersion flags
    return old_edisp;
}


/***********************************************************************//**
 * @brief Restore energy dispersion flags of CTA observations
 *
 * @param[in] obs Observation container.
 * @param[in] edisp Vector of energy dispersion flags.
 *
 * Restores energy dispersion flags that have previously been returned by
 * the set_edisp() method. The number of observations between the calls
 * of the set_edisp() and this method has to be the same.
 ***************************************************************************/
void ctool::restore_edisp(GObservations& obs, const std::vector<bool>& edisp) const
{
    // Check that energy dispersion flag vector has the same size as the
    // observation container
    if (edisp.size() != obs.size()) {
        std::string msg = "Size of energy dispersion flag vector ("+
                          gammalib::str(edisp.size())+" elements) is "
                          "incompatible with number of observations in the "
                          "observation container ("+gammalib::str(obs.size())+
                          " observation).";
        throw GException::invalid_value(G_RESTORE_EDISP, msg);
    }

    // Loop over all observations in observation container
    for (int i = 0; i < obs.size(); ++i) {

        // Is this observation a CTA observation?
        GCTAObservation* cta = dynamic_cast<GCTAObservation*>(obs[i]);

        // Yes, then restore the energy dispersion flag
        if (cta != NULL) {
            cta->response()->apply_edisp(edisp[i]);
        }

    } // endfor: loop over all observations in container

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set response for CTA observation
 *
 * @param[in,out] obs CTA observation
 *
 * @exception GException::invalid_value
 *            Energy dispersion requested but not energy dispersion cube given
 *
 * Set the response for one CTA observation. If the CTA observation contains
 * a counts cube the method first attempts to set the response information
 * using the following user parameters
 *
 *      expcube:   Exposure cube file
 *      psfcube:   Point spread function cube file
 *      edispcube: Energy dispersion cube file (optional)
 *      bkgcube:   Background cube file
 *
 * If none of the @p expcube, @p psfcube and @p bkgcube parameters is
 * "NONE" or empty, the method will allocate a stacked response for the
 * observation. In case that an energy dispersion cube is provided using
 * the @p edispcube parameter, this energy dispersion cube is also loaded
 * and attached to the stacked response function.
 *
 * In case that any of the @p expcube, @p psfcube and @p bkgcube parameters
 * is "NONE" or empty, the method will assume that a IRF response should
 * be used, and will attempty to set the response information using the
 * following user parameters
 *
 *      caldb:     Calibration database
 *      irf:       Instrument response function
 *
 ***************************************************************************/
void ctool::set_obs_response(GCTAObservation* obs)
{
    // Initialise response flag
    bool has_response = false;

    // If the observation contains a counts cube, then check whether
    // response information for a stacked analysis is provided, and if so,
    // set that response information for the CTA observation
    if (dynamic_cast<const GCTAEventCube*>(obs->events()) != NULL) {

        // Check if response cube parameters are provided. We need all four
        // response cube parameters.
        if (has_par("expcube") && has_par("psfcube") &&
            has_par("bkgcube") && has_par("edispcube")) {

            // Check if we need to query the energy dispersion cube. We do
            // not need to query the cube is an "edisp" parameter exists and
            // if this parameter is set to "no".
            bool query_edisp = true;
            if (has_par("edisp")) {
                query_edisp = (*this)["edisp"].boolean();
            }

            // Get filenames
            GFilename expcube((*this)["expcube"].filename());
            GFilename psfcube((*this)["psfcube"].filename());
            GFilename edispcube;
            if (query_edisp) {
                edispcube = (*this)["edispcube"].filename();
            }
            GFilename bkgcube((*this)["bkgcube"].filename());

            // Extract stacked response information if available
            if (is_valid_filename(expcube) && is_valid_filename(psfcube) &&
                is_valid_filename(bkgcube)) {

                // Load exposure, PSF and background cubes
                GCTACubeExposure   exposure(expcube);
                GCTACubePsf        psf(psfcube);
                GCTACubeBackground background(bkgcube);

                // If querying of energy dispersion cube is requested then
                // query it now
                if (query_edisp) {

                    // If filename is valid then use energy dispersion ...
                    if (is_valid_filename(edispcube)) {

                        // Load energy dispersion cube
                        GCTACubeEdisp edisp(edispcube);

                        // Set response with all four cubes
                        obs->response(exposure, psf, edisp, background);

                    } // endif: energy dispersion cube was provided

                    // ... otherwise work without energy dispersion
                    else {

                        // If energy dispersion was requested but no cube
                        // was provided then throw an exception
                        if (has_par("edisp") && (*this)["edisp"].boolean()) {
                            std::string msg = "Energy dispersion requested but "
                                              "no energy dispersion cube was "
                                              "specified.";
                            throw GException::invalid_value(G_SET_OBS_RESPONSE,
                                                            msg);
                        }

                        // Set response without energy dispersion
                        obs->response(exposure, psf, background);
                        
                    } // endelse: no energy dispersion was available

                } // endif: energy dispersion needed

                // ... otherwise work without energy dispersion
                else {
                    obs->response(exposure, psf, background);
                }

                // Signal that response is available
                has_response = true;

            } // endif: filenames were available

        } // endif: expcube, psfcube and bkdcube parameters valid
    } // endif: observation contains a counts cube

    // If we have not yet response information then get it now from the
    // caldb and irf user parameters
    if (!has_response) {

        // Load response information
        std::string database = (*this)["caldb"].string();
        std::string irf      = (*this)["irf"].string();

        // Create an XML element containing the database and IRF name. This
        // kluge will make sure that the information is later written
        // to the observation definition XML file, in case an observation
        // definition XML file is written.
        std::string parameter = "parameter name=\"Calibration\""
                                " database=\""+database+"\""
                                " response=\""+irf+"\"";

        GXmlElement xml;
        xml.append(parameter);

        // Create CTA response
        GCTAResponseIrf response(xml);

        // Attach response to observation
        obs->response(response);
        
        // Signal response availability
        has_response = true;

    } // endif: no response information was available

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set observation boundaries for all CTA observations
 *
 * @param[in,out] obs Observation container
 *
 * Sets the observation boundaries for all CTA observations in the container
 * that contain event lists. If the energy boundaries of an event list are
 * empty, the method gets the 
 *
 *      emin:   Minimum energy (TeV)
 *      emax:   Maximum energy (TeV)
 *
 * user parameters and sets a single energy boundary element @p [emin,emax].
 *
 * If the RoI radius of the event list is invalid, the method gets the
 *
 *      rad:    Region of Interest radius (deg)
 *
 * user parameter to set the radius of the RoI, and uses the CTA pointing
 * direction to set the centre of the RoI.
 ***************************************************************************/
void ctool::set_obs_bounds(GObservations& obs)
{
    // Loop over all observations in the container
    for (int i = 0; i < obs.size(); ++i) {

        // Get pointer to CTA observation
        GCTAObservation* cta = dynamic_cast<GCTAObservation*>(obs[i]);

        // Fall through if observation is not a CTA observation or if it
        // does not contain events
        if ((cta == NULL) || (!cta->has_events())) {
            continue;
        }

        // Get pointer on CTA event list
        GCTAEventList* list = const_cast<GCTAEventList*>
                       (dynamic_cast<const GCTAEventList*>(cta->events()));

        // Fall through if CTA observation does not contain an event list
        if (list == NULL) {
            continue;
        }

        // If there are no energy boundaries then read the "emin" and "emax"
        // user parameters and add them
        if (list->ebounds().is_empty()) {
            if (has_par("emin") && has_par("emax")) {
                double emin((*this)["emin"].real());
                double emax((*this)["emax"].real());
                GEbounds ebounds(GEnergy(emin, "TeV"),
                                 GEnergy(emax, "TeV"));
                list->ebounds(ebounds);
            }
        }

        // If there is no RoI then read the "rad" user parameters and use
        // the pointing direction to set the RoI
        if (list->roi().radius() == 0) {
            if (has_par("rad")) {
                double rad = (*this)["rad"].real();
                GCTARoi roi(GCTAInstDir(cta->pointing().dir()), rad);
                list->roi(roi);
            }
        }
        
    } // endfor: looped over observations

    // Return
    return;
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
    if (!is_valid_filename(filename)) {

        // Setup a new CTA observation
        GCTAObservation cta_obs = create_cta_obs();

        // Get response if required
        if (get_response) {

            // Set response
            set_obs_response(&cta_obs);

        } // endif: response was required

       // Append observation to container
       obs.append(cta_obs);

    } // endif: filename was valid

    // ... otherwise we have a file name
    else {

        // If file is a FITS file then create an empty CTA observation
        // and load file into observation
        if (GFilename(filename).is_fits()) {

            // Allocate empty CTA observation
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

        // ... otherwise load file into observation container
        else {

            // Load observations from XML file
            obs.load(filename);

            // Get response if required
            if (get_response) {

                // For all observations that have no response, set the response
                // from the task parameters
                set_response(obs);

            } // endif: response was required

            // Set observation boundary parameters (emin, emax, rad)
            set_obs_bounds(obs);

            // Signal that XML file should be used for storage
            m_use_xml = true;

        } // endelse: file was an XML file

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
 * @brief Get current resident set size (physical memory use) in Bytes
 *
 * @return Physical memory use in Bytes.
 *
 * @todo This method is currently not used and can be removed if we do not
 *       use it
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
 * @brief Return observation header string
 *
 * @param[in] obs Pointer to observation.
 * @return String containing observation information.
 *
 * Returns a string that contains observation information, including the
 * instrument name, the observation name, and the observation ID. The format
 * of the string is
 *
 *      XXX observation "name" (id=YYYYY)
 *
 ***************************************************************************/
std::string ctool::get_obs_header(const GObservation* obs) const
{
    // Initialise header string
    std::string header = obs->instrument() + " observation";

    // If observation name is not empty then add name
    if (!obs->name().empty()) {
        header += " \"" + obs->name() + "\"";
    }

    // If observation ID is not empty then add ID
    if (!obs->id().empty()) {
        header += " (id=" + obs->id() +")";
    }

    // Return header
    return header;
}


/***********************************************************************//**
 * @brief Insert observation energy boundaries into list of energies
 *
 * @param[in] energies Energies.
 * @param[in] obs Observation container.
 * @return Energies.
 *
 * Inserts the energy boundaries of an observation in a list of @p energies.
 ***************************************************************************/
GEnergies ctool::insert_energy_boundaries(const GEnergies&       energies,
                                          const GCTAObservation& obs)
{
    // Create copy of input energies
    GEnergies engs = energies;

    // Get the energy boundaries of the event list
    GEbounds ebounds = obs.ebounds();

    // Loop over all boundary energies
    for (int iebin = 0; iebin <= ebounds.size(); ++iebin) {

        // Get boundary energy
        GEnergy energy;
        if (iebin < ebounds.size()) {
            energy = ebounds.emin(iebin);
        }
        else {
            energy = ebounds.emax();
        }

        // Get index before which the energy should be appended
        bool insert = true;
        int  index  = -1;
        for (int k = 0; k < engs.size(); ++k) {
        
            // If energy exists already then skip the boundary and examine the
            // next one. We consider here 1 MeV as being sufficiently close.
            if (std::abs(engs[k].MeV()-energy.MeV()) < 1.0) {
                insert = false;
                break;
            }
        
            // If energy is above the boundary energy we found the index
            if (engs[k] > energy) {
                index = k;
                break;
            }

        } // endfor: search index for insertion

        // Insert energy if requested
        if (insert) {

            // Insert or append energy
            if (index != -1) {
                engs.insert(index, energy);
            }
            else {
                engs.append(energy);
            }

            // Log energy insertion. Circumvent const correctness
            log_value(NORMAL, "Insert energy", energy.print());

        } // endif: energy insertion requested

    } // endfor: looped over all energy boundaries

    // Return energies
    return engs;
}


/***********************************************************************//**
 * @brief Determine the counts cube layer usage.
 *
 * @param[in] cube_ebounds Energy boundaries of the counts cube.
 * @param[in] list_ebounds Energy boundaries of the event list.
 * @return Vector of usage flags.
 *
 * Determines a vector of counts cube layer usage flags that signal whether
 * an event should actually be filled in the counts cube or not. This makes
 * sure that no partially filled bins will exist in the counts cube.
 ***************************************************************************/
std::vector<bool> ctool::cube_layer_usage(const GEbounds& cube_ebounds,
                                          const GEbounds& list_ebounds) const
{
    // Set energy margin
    const GEnergy energy_margin(0.01, "GeV");

    // Initialise usage vector
    std::vector<bool> usage(cube_ebounds.size(), true);

    // Loop over all energy bins of the cube
    for (int iebin = 0; iebin < cube_ebounds.size(); ++iebin) {

        // If the counts cube energy bin is fully contained within the
        // energy boundaries of the event list the signal usage of this
        // counts cube bin. Partially overlapping energy bins are signaled
        // for non-usage, which avoids having partially filled bins in the
        // counts cube. Some margin is applied that effectively reduces the
        // width of the counts cube energy bin. This should cope with any
        // rounding errors that come from reading and writing the energy
        // boundary information to a FITS file.
        usage[iebin] = list_ebounds.contains(cube_ebounds.emin(iebin) +
                                             energy_margin,
                                             cube_ebounds.emax(iebin) -
                                             energy_margin);

    } // endfor: looped over all energy bins

    // Return usage
    return usage;
}


/***********************************************************************//**
 * @brief Set warning string if there are too few energies
 *
 * @param[in] energies Energies.
 * @return Warning string.
 *
 * Sets a warning string if there are too few @p energies in the container.
 * For a binned or stacked analysis to be accurate, at least 25 bins per
 * decade are required. If the provided number of energies is less than this
 * number a warning string will be returned. Otherwise an empty string will
 * be returned.
 ***************************************************************************/
std::string ctool::warn_too_few_energies(const GEnergies& energies) const
{
    // Initialise warning
    std::string warning;

    // Retrieve the number of energies
    int n = energies.size();

    // Compute the required number of energy bins
    double logEmin   = std::log10(energies[0].TeV());
    double logEmax   = std::log10(energies[n-1].TeV());
    int    nrequired = int((logEmax - logEmin) * 25.0);
    
    // Warn if there are not enough energy bins
    if (n < nrequired) {
        warning.append("\nWARNING: Only "+gammalib::str(n-1)+" energy bins "
                       "have been requested. This may be too few energy "
                       "bins."
                       "\n         At least 25 bins per decade in energy are "
                       "recommended, which for the given"
                       "\n         energy range would be "+
                       gammalib::str(nrequired)+" bins. Consider increasing "
                       "the number of energy bins.");
    }

    // Return warning
    return warning;
}


/***********************************************************************//**
 * @brief Set warning string if file has no .xml suffix
 *
 * @param[in] filename Filename.
 * @return Warning string.
 *
 * Sets a warning string if the filename has not ".xml" suffix.
 ***************************************************************************/
std::string ctool::warn_xml_suffix(const GFilename& filename) const
{
    // Initialise warning
    std::string warning;

    // Get filename suffix
    std::string fname  = std::string(filename);
    std::string suffix = gammalib::tolower(fname.substr(fname.length()-4,4));
    
    // Warn if there are not enough energy bins
    if (suffix != ".xml") {
        warning.append("WARNING: Name of observation definition output file "
                       "\""+filename+"\" does not terminate with \".xml\"."
                       "\n         This is not an error, but may be misleading. "
                       "It is recommended to use the suffix"
                       "\n         \".xml\" for observation definition files.");
    }

    // Return warning
    return warning;
}


/*==========================================================================
 =                                                                         =
 =                         Private methods for SWIG                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Dump help text in the console
 *
 * Dumps the help text for the ctool into the console. The help text is
 * located in the folder
 *
 *      $CTOOLS//share/help/
 *
 * of the ctools installation and the help file has the name "name().txt",
 * where name() stands for the name of the ctool.
 ***************************************************************************/
void ctool::provide_help(void) const
{
    // Allocate line buffer
    const int n = 1000; 
    char  line[n];

    // Build help file name
    std::string helpfile = name()+".txt";

    // Get ctools environment variable
    char* ptr = std::getenv("CTOOLS");
    if (ptr == NULL) {
        std::string msg = "CTOOLS environment variable not set, cannot "
                          "display help file. Please set the CTOOLS "
                          "environment variable.";
        throw GException::invalid_value(G_PROVIDE_HELP, msg);
    }

    // If help file exists then display it, otherwise notify that no
    // help is available
    std::string fname = std::string(ptr) + "/share/help/" + helpfile;
    FILE* fptr = fopen(fname.c_str(), "r");
    if (fptr != NULL) {
        while (fgets(line, n, fptr) != NULL) {
            std::cout << std::string(line);
        }
        fclose(fptr);
    }
    else {
        std::cout << "No help available for "+name()+"." << std::endl;
    }

    // Return
    return;
}
