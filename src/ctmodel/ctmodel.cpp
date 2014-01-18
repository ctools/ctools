/***************************************************************************
 *                     ctmodel - CTA counts model tool                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @file ctmodel.cpp
 * @brief CTA counts model tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctmodel.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_SETUP_OBS                                    "ctmodel::setup_obs()"
#define G_MODEL_MAP                    "ctmodel::model_map(GCTAObservation*)"

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
ctmodel::ctmodel(void) : GApplication(CTMODEL_NAME, CTMODEL_VERSION)
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
ctmodel::ctmodel(GObservations obs) : GApplication(CTMODEL_NAME, CTMODEL_VERSION)
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
ctmodel::ctmodel(int argc, char *argv[]) :
         GApplication(CTMODEL_NAME, CTMODEL_VERSION, argc, argv)
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
ctmodel::ctmodel(const ctmodel& app) : GApplication(app)
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
ctmodel::~ctmodel(void)
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
 * @param[in] app ctmodel application.
 * @return Returns ctmodel application.
 ***************************************************************************/
ctmodel& ctmodel::operator= (const ctmodel& app)
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
void ctmodel::clear(void)
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
 * This is the main execution method of the ctmodel class. It is invoked
 * when the executable is called from command line. The method generates
 * the model maps and saves the results.
 ***************************************************************************/
void ctmodel::execute(void)
{
    // Signal that some parameters should be read ahead
    m_read_ahead = true;

    // Create the model map(s)
    run();

    // Save the model map(s) into FITS file
    save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generate the model map(s)
 *
 * This method reads the task parameters from the parfile, sets up the
 * observation container, loops over all CTA observations in the container
 * and generates a model map for each CTA observation.
 ***************************************************************************/
void ctmodel::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Setup observation container
    setup_obs();

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
            log.header1("Generate model maps");
        }
        else {
            log.header1("Generate model map");
        }
    }

    // Initialise observation counter
    int n_observations = 0;

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

            // Generate model map
            model_map(obs, m_obs.models());

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
            log.header1("Observations after model map generation");
        }
        else {
            log.header1("Observation after model map generation");
        }
        log << m_obs << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save model map(s)
 *
 * This method saves the model map(s) into FITS file(s). There are two
 * modes, depending on the m_use_xml flag.
 *
 * If m_use_xml is true, all model map(s) will be saved into FITS files,
 * where the output filenames are constructued from the input filenames by
 * prepending the m_prefix string to name. Any path information will be
 * stripped form the input name, hence event files will be written into the
 * local working directory (unless some path information is present
 * in the prefix). In addition, an XML file will be created that gathers
 * the filename information for the model map(s). If an XML file was present
 * on input, all metadata information will be copied from this input file.
 *
 * If m_use_xml is false, the model map will be saved into a FITS file.
 ***************************************************************************/
void ctmodel::save(void)
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

    // Case A: Save model map(s) and XML metadata information
    if (m_use_xml) {
        save_xml();
    }

    // Case B: Save model map as FITS file
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
 * the user. The parameters are read in the correct order.
 ***************************************************************************/
void ctmodel::get_parameters(void)
{
    // If we do not have any observations in the container then get all
    // information to build (at least) one
    if (m_obs.size() == 0) {

        // Read general parameters. Only read the output parameters if
        // required, and also only read the model XML filename if no
        // model is yet in the container.
        m_infile = (*this)["infile"].filename();
        if (m_read_ahead) {
            m_outfile = (*this)["outfile"].filename();
            m_prefix  = (*this)["prefix"].string();
        }
        m_caldb  = (*this)["caldb"].string();
        m_irf    = (*this)["irf"].string();
        if (m_obs.models().size() == 0) {
            m_srcmdl = (*this)["srcmdl"].filename();
        }

        // If there is no input filename then read all parameters that
        // are required to build a model map from scratch
        if ((m_infile == "NONE") || (gammalib::strip_whitespace(m_infile) == "")) {
            m_ra       = (*this)["ra"].real();
            m_dec      = (*this)["dec"].real();
            m_deadc    = (*this)["deadc"].real();
            m_tmin     = (*this)["tmin"].real();
            m_tmax     = (*this)["tmax"].real();
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
        }
    }

    // ... otherwise, read only the parameters that are required for the
    // model generation.
    else {
        if (m_read_ahead) {
            m_outfile = (*this)["outfile"].filename();
            m_prefix  = (*this)["prefix"].string();
        }
        m_caldb  = (*this)["caldb"].string();
        m_irf    = (*this)["irf"].string();
        if (m_obs.models().size() == 0) {
            m_srcmdl = (*this)["srcmdl"].filename();
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup observation container
 *
 * @exception GException::no_cube
 *            No event cube found in CTA observation.
 *
 * This method sets up the observation container for processing. There are
 * two cases:
 *
 * If there are no observations in the actual observation container, the
 * method will check in "infile" parameter. If this parameter is "NONE" or
 * empty, the task parameters will be used to construct a model map.
 * Otherwise, the method first tries to interpret the "infile" parameter as
 * a counts map, and attemps loading of the file in an event cube. If this
 * fails, the method tries to interpret the "infile" parameter as an
 * observation definition XML file. If this also fails, an exception will
 * be thrown.
 *
 * If observations exist already in the observation container, the method
 * will simply keep them.
 *
 * Test if all CTA observations contain counts maps.
 *
 * Finally, if no models exist so far in the observation container, the
 * models will be loaded from the model XML file.
 ***************************************************************************/
void ctmodel::setup_obs(void)
{
    // If there are no observations in the container then try to build some
    if (m_obs.size() == 0) {
        
        // If no input filename has been specified, then create a model map
        // from the task parameters
        if ((m_infile == "NONE") || (gammalib::strip_whitespace(m_infile) == "")) {

            // Set pointing direction
            GCTAPointing pnt;
            GSkyDir      skydir;
            skydir.radec_deg(m_ra, m_dec);
            pnt.dir(skydir);

            // Setup energy range covered by model
            GEnergy  emin(m_emin, "TeV");
            GEnergy  emax(m_emax, "TeV");
            GEbounds ebds(m_enumbins, emin, emax);

            // Setup time interval covered by model
            GGti  gti;
            GTime tmin(m_tmin);
            GTime tmax(m_tmax);
            gti.append(tmin, tmax);

            // Setup skymap
            GSkymap map = GSkymap(m_proj, m_coordsys,
                                  m_xref, m_yref, -m_binsz, m_binsz,
                                  m_nxpix, m_nypix, m_enumbins);

            // Create model cube from sky map
            GCTAEventCube cube(map, ebds, gti);

            // Allocate CTA observation
            GCTAObservation obs;

            // Set CTA observation attributes
            obs.pointing(pnt);
            obs.ontime(gti.ontime());
            obs.livetime(gti.ontime()*m_deadc);
            obs.deadc(m_deadc);

            // Set event cube in observation
            obs.events(cube);

            // Append CTA observation to container
            m_obs.append(obs);

            // Signal that no XML file should be used for storage
            m_use_xml = false;

        } // endif: created model map from task parameters

        // ... otherwise try to load information from the file
        else {

            // First try to open the file as a counts map
            try {

                // Allocate CTA observation
                GCTAObservation obs;

                // Load counts map in CTA observation
                obs.load_binned(m_infile);

                // Append CTA observation to container
                m_obs.append(obs);

                // Signal that no XML file should be used for storage
                m_use_xml = false;
            
            }
        
            // ... otherwise try to open as XML file
            catch (GException::fits_open_error &e) {

                // Load observations from XML file. This will throw
                // an exception if it fails.
                m_obs.load(m_infile);

                // Signal that XML file should be used for storage
                m_use_xml = true;

            }

        } // endelse: loaded information from input file

    } // endif: there was no observation in the container

    // If there are no models associated with the observations then
    // load now the model definition from the XML file
    if (m_obs.models().size() == 0) {
        m_obs.models(GModels(m_srcmdl));
    }

    // Check if all CTA observations contain an event cube and setup response
    // for all observations
    for (int i = 0; i < m_obs.size(); ++i) {

        // Is this observation a CTA observation?
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Yes ...
        if (obs != NULL) {

            // Throw an exception if this observation does not contain
            // an event cube
            if (dynamic_cast<const GCTAEventCube*>(obs->events()) == NULL) {
                throw GException::no_cube(G_SETUP_OBS);
            }
        
            // Set response if it isn't set already
            if (obs->response().aeff() == NULL) {
            	obs->response(m_irf, m_caldb);

            } // endif: observation already has a response

        } // endif: observation was a CTA observation

    } // endfor: looped over all observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generate model map
 *
 * @param[in] obs CTA observation pointer.
 * @param[in] models Model container.
 *
 * @exception GException::no_cube
 *            No event cube found in CTA observation.
 ***************************************************************************/
void ctmodel::model_map(GCTAObservation* obs, const GModels& models)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Get event cube pointer
        GCTAEventCube* cube = 
            const_cast<GCTAEventCube*>(dynamic_cast<const GCTAEventCube*>(obs->events()));

        // Throw an exception if the observation does not hold and event
        // cube
        if (cube == NULL) {
            throw GException::no_cube(G_MODEL_MAP);
        }

        // Initialise statistics
        double sum = 0.0;

        // Loop over all events in counts map
        for (int i = 0; i < cube->size(); ++i) {

            // Get event bin
            GCTAEventBin* bin = (*cube)[i];
            
            // Compute model value for event bin
            double model = 
                   models.eval(*(const_cast<const GCTAEventBin*>(bin)), *obs) *
                   bin->size();

            // Store value
            bin->counts(model);

            // Sum all events
            sum += model;
        }

        // Log results
        if (logTerse()) {
            log << gammalib::parformat("Model events in cube");
            log << sum << std::endl;
        }

    } // endif: observation pointer was not valid

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
void ctmodel::init_members(void)
{
    // Initialise members
    m_infile.clear();
    m_outfile.clear();
    m_prefix.clear();
    m_caldb.clear();
    m_irf.clear();
    m_srcmdl.clear();
    m_proj.clear();
    m_coordsys.clear();
    m_ra       = 0.0;
    m_dec      = 0.0;
    m_deadc    = 1.0;
    m_tmin     = 0.0;
    m_tmax     = 0.0;
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
void ctmodel::copy_members(const ctmodel& app)
{
    // Copy attributes
    m_infile   = app.m_infile;
    m_outfile  = app.m_outfile;
    m_prefix   = app.m_prefix;
    m_caldb    = app.m_caldb;
    m_irf      = app.m_irf;
    m_srcmdl   = app.m_srcmdl;
    m_proj     = app.m_proj;
    m_coordsys = app.m_coordsys;
    m_ra       = app.m_ra;
    m_dec      = app.m_dec;
    m_deadc    = app.m_deadc;
    m_tmin     = app.m_tmin;
    m_tmax     = app.m_tmax;
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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctmodel::free_members(void)
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
 * Converts an input filename into an output filename by prepending the
 * prefix stored in the member m_prefix to the input filename. Any path will
 * be stripped from the input filename. Also a trailing ".gz" will be
 * stripped.
 ***************************************************************************/
std::string ctmodel::set_outfile_name(const std::string& filename) const
{
    // Split input filename into path elements
    std::vector<std::string> elements = gammalib::split(filename, "/");

    // The last path element is the filename
    std::string outname = m_prefix + elements[elements.size()-1];

    // Strip any ".gz"
    outname = gammalib::strip_chars(outname, ".gz");
    
    // Return output filename
    return outname;
}


/***********************************************************************//**
 * @brief Save model map in FITS format.
 *
 * Save the model map as a FITS file. The filename of the FITS file is
 * specified by the m_outfile member.
 ***************************************************************************/
void ctmodel::save_fits(void)
{
    // Get output filename
    m_outfile = (*this)["outfile"].filename();

    // Get CTA observation from observation container
    GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[0]);

    // Save model map
    save_model_map(obs, m_outfile);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save model map(s) in XML format.
 *
 * Save the model map(s) into FITS files and write the file path information
 * into a XML file. The filename of the XML file is specified by the
 * m_outfile member, the filename(s) of the model map(s) are built by
 * prepending the prefix given by the m_prefix member to the input model
 * map(s) filenames. Any path present in the input filename will be stripped,
 * i.e. the model map(s) will be written in the local working directory
 * (unless a path is specified in the m_prefix member).
 ***************************************************************************/
void ctmodel::save_xml(void)
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
            std::string outfile = set_outfile_name(m_infiles[i]);

            // Store output file name in observation
            obs->eventfile(outfile);

            // Save event list
            save_model_map(obs, outfile);

        } // endif: observation was a CTA observations

    } // endfor: looped over observations

    // Save observations in XML file
    m_obs.save(m_outfile);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save a single model map into a FITS file
 *
 * @param[in] obs Pointer to CTA observation.
 * @param[in] outfile Output file name.
 *
 * This method saves a single model map into a FITS file. The method does
 * nothing if the observation pointer is not valid.
 ***************************************************************************/
void ctmodel::save_model_map(const GCTAObservation* obs,
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
