/***************************************************************************
 *                  ctmodel - Model cube generation tool                   *
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
 * @brief Model cube generation tool implementation
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
ctmodel::ctmodel(void) : ctool(CTMODEL_NAME, CTMODEL_VERSION)
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
ctmodel::ctmodel(const GObservations& obs) : ctool(CTMODEL_NAME, CTMODEL_VERSION)
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
ctmodel::ctmodel(int argc, char *argv[]) :
         ctool(CTMODEL_NAME, CTMODEL_VERSION, argc, argv)
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
ctmodel::ctmodel(const ctmodel& app) : ctool(app)
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
 * @param[in] app Application.
 * @return Application.
 ***************************************************************************/
ctmodel& ctmodel::operator=(const ctmodel& app)
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
void ctmodel::clear(void)
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

    // Set energy dispersion flag for all CTA observations and save old
    // values in save_edisp vector
    std::vector<bool> save_edisp;
    save_edisp.assign(m_obs.size(), false);
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            save_edisp[i] = obs->response()->apply_edisp();
            obs->response()->apply_edisp(m_apply_edisp);
        }
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
        log.header1("Generate model cube");
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

    // Log cube
    if (logTerse()) {
        log << std::endl;
        log.header1("Model cube");
        log << m_cube << std::endl;
    }

    // Restore energy dispersion flag for all CTA observations
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            obs->response()->apply_edisp(save_edisp[i]);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save model cube
 *
 * Saves the model cube into a FITS file specified using the "outfile"
 * task parameter.
 ***************************************************************************/
void ctmodel::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save cube");
    }

    // Make sure we have the FITS filename
    m_outfile = (*this)["outfile"].filename();

    // Save model cube into FITS file
    m_cube.save(m_outfile, clobber());

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
    m_obsfile.clear();
    m_outfile.clear();
    m_caldb.clear();
    m_irf.clear();
    m_srcmdl.clear();
    m_proj.clear();
    m_coordsys.clear();
    m_ebinalg.clear();
    m_ebinfile.clear();
    m_ra          = 0.0;
    m_dec         = 0.0;
    m_deadc       = 1.0;
    m_tmin        = 0.0;
    m_tmax        = 0.0;
    m_emin        = 0.0;
    m_emax        = 0.0;
    m_enumbins    = 0;
    m_xref        = 0.0;
    m_yref        = 0.0;
    m_binsz       = 0.0;
    m_nxpix       = 0;
    m_nypix       = 0;
    m_apply_edisp = false;

    // Initialise protected members
    m_obs.clear();
    m_cube.clear();
    m_ebounds.clear();
    m_gti.clear();

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
    m_infile      = app.m_infile;
    m_obsfile     = app.m_obsfile;
    m_outfile     = app.m_outfile;
    m_caldb       = app.m_caldb;
    m_irf         = app.m_irf;
    m_srcmdl      = app.m_srcmdl;
    m_proj        = app.m_proj;
    m_coordsys    = app.m_coordsys;
    m_ebinalg     = app.m_ebinalg;
    m_ebinfile    = app.m_ebinfile;
    m_ra          = app.m_ra;
    m_dec         = app.m_dec;
    m_deadc       = app.m_deadc;
    m_tmin        = app.m_tmin;
    m_tmax        = app.m_tmax;
    m_emin        = app.m_emin;
    m_emax        = app.m_emax;
    m_enumbins    = app.m_enumbins;
    m_xref        = app.m_xref;
    m_yref        = app.m_yref;
    m_binsz       = app.m_binsz;
    m_nxpix       = app.m_nxpix;
    m_nypix       = app.m_nypix;
    m_apply_edisp = app.m_apply_edisp;

    // Copy protected members
    m_obs        = app.m_obs;
    m_cube       = app.m_cube;
    m_ebounds    = app.m_ebounds;
    m_gti        = app.m_gti;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctmodel::free_members(void)
{
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
    // Read input and optionally output cube filenames
    m_infile = (*this)["infile"].filename();
    if (read_ahead()) {
        m_outfile = (*this)["outfile"].filename();
    }

    // Read model definiton file if required
    if (m_obs.models().size() == 0) {
        m_srcmdl = (*this)["srcmdl"].filename();
    }

    // Get energy dispersion flag parameters
    m_apply_edisp = (*this)["edisp"].boolean();

    // If the observation container is empty then read the parameters that
    // are necessary to fill it
    if (m_obs.size() == 0) {
    
        // Read observation definition filename and response parameters
        m_obsfile = (*this)["obsfile"].filename();
        m_caldb   = (*this)["caldb"].string();
        m_irf     = (*this)["irf"].string();

        // If no observation definition file has been specified then read all
        // parameters that are necessary to create an observation from scratch
        // (see method setup_obs)
        if ((m_obsfile == "NONE") || (gammalib::strip_whitespace(m_obsfile) == "")) {
            m_ra    = (*this)["ra"].real();
            m_dec   = (*this)["dec"].real();
            m_deadc = (*this)["deadc"].real();
            m_tmin  = (*this)["tmin"].real();
            m_tmax  = (*this)["tmax"].real();
            m_emin  = (*this)["emin"].real();
            m_emax  = (*this)["emax"].real();
        }

    }

    // Check if we need response information
    for (int i = 0; i < m_obs.size(); ++i) {
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);
        if (obs != NULL) {
            if (!obs->has_response()) {
                m_caldb = (*this)["caldb"].string();
                m_irf   = (*this)["irf"].string();
                break;
            }
        }
    }

    // If no cube file has been specified then read all parameters that
    // are necessary to create the cube from scratch (see method init_cube)
    if ((m_infile == "NONE") || (gammalib::strip_whitespace(m_infile) == "")) {
        m_ebinalg  = (*this)["ebinalg"].string();
        if (m_ebinalg == "FILE") {
            m_ebinfile = (*this)["ebinfile"].filename();
        }
        else {
            m_emin     = (*this)["emin"].real();
            m_emax     = (*this)["emax"].real();
            m_enumbins = (*this)["enumbins"].integer();
        }
        m_proj     = (*this)["proj"].string();
        m_coordsys = (*this)["coordsys"].string();
        m_xref     = (*this)["xref"].real();
        m_yref     = (*this)["yref"].real();
        m_binsz    = (*this)["binsz"].real();
        m_nxpix    = (*this)["nxpix"].integer();
        m_nypix    = (*this)["nypix"].integer();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup observation container
 *
 * This method makes sure that we have an observation container with an
 * associated model and that all CTA observations have a valid response
 * functions.
 *
 * If on input the observation container is empty, the method sets up an
 * observation container by either loading the information from an
 * observation definition XML file, by loading an event list or a counts
 * cube as a single observation, or by constructing an observation from
 * the scratch from the task parameters.
 *
 * If no model exists yet in the observation container, the method will load
 * the model from the file given by the task parameter "srcmdl".
 *
 * The method then loops over all observations in the container and makes
 * sure that all CTA observations have a valid response function. If no
 * response function exists so far, the response function will be set
 * using the "caldb" and "irf" task parameters.
 ***************************************************************************/
void ctmodel::setup_obs(void)
{
    // If there are no observations in the container then get them either
    // from the observation definition file or from the task parameters
    if (m_obs.size() == 0) {
        
        // If no observation definition file has been specified, then create
        // an observation container from the task parameters
        if ((m_obsfile == "NONE") || (gammalib::strip_whitespace(m_obsfile) == "")) {

            // Set pointing direction
            GCTAPointing pnt;
            GSkyDir      skydir;
            skydir.radec_deg(m_ra, m_dec);
            pnt.dir(skydir);

            // Setup time interval covered by model
            GGti  gti;
            GTime tmin(m_tmin);
            GTime tmax(m_tmax);
            gti.append(tmin, tmax);

            // Setup energy range covered by model
            GEnergy  emin(m_emin, "TeV");
            GEnergy  emax(m_emax, "TeV");
            GEbounds ebounds(1, emin, emax);

            // Allocate CTA observation and empty event list
            GCTAObservation obs;
            GCTAEventList   list;

            // Set event list GTI and energy boundaries and append event list
            // to observation (we need in fact only the GTI for the
            // computations)
            list.gti(gti);
            list.ebounds(ebounds);
            obs.events(list);

            // Set CTA observation attributes
            obs.pointing(pnt);
            obs.ontime(gti.ontime());
            obs.livetime(gti.ontime()*m_deadc);
            obs.deadc(m_deadc);

            // Append CTA observation to container
            m_obs.append(obs);

        } // endif: created single CTA observation from task parameters

        // ... otherwise try to load information from the file
        else {

            // First try to open the file as an event list or counts map
            try {

                // Allocate CTA observation
                GCTAObservation obs;

                // Load CTA observation from FITS file
                obs.load(m_obsfile);

                // Append CTA observation to container
                m_obs.append(obs);
            
            }
        
            // ... otherwise try to open as XML file
            catch (GException::fits_open_error &e) {

                // Load observations from XML file. This will throw
                // an exception if it fails.
                m_obs.load(m_obsfile);

            }

        } // endelse: loaded information from input file

    } // endif: there was no observation in the container

    // If there are no models associated with the observations then
    // load now the model definition from the XML file
    if (m_obs.models().size() == 0) {
        m_obs.models(GModels(m_srcmdl));
    }

    // For all observations that have no response, set the response
    // from the task parameters
    set_response(m_obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise model cube
 *
 * Initialises model cube by either reading the cube definition from an
 * existing cube file specified by the "infile" parameter, or by constructing
 * a cube from the energy and spatial binning information provided by the
 * task parameters. The GTIs of the model cube will be undefined on
 * exit.
 ***************************************************************************/
void ctmodel::init_cube(void)
{
    // Clear model cube and GTIs
    m_cube.clear();
    m_gti.clear();

    // If no cube file has been specified then create a cube from the task
    // parameters
    if ((m_infile == "NONE") || (gammalib::strip_whitespace(m_infile) == "")) {

        // Set dummy GTI that is needed for event cube creation
        GGti gti;
        gti.append(GTime(0.0), GTime(0.1234));
    
        // Set energy boundaries
        m_ebounds = get_ebounds();

        // Setup skymap
        GSkymap map = GSkymap(m_proj, m_coordsys,
                              m_xref, m_yref, -m_binsz, m_binsz,
                              m_nxpix, m_nypix, m_enumbins);

        // Create model cube
        m_cube = GCTAEventCube(map, m_ebounds, gti);

    }

    // ... otherwise load cube from file
    else {

        // Load cube
        m_cube.load(m_infile);

        // Set all cube bins to zero
        for (int i = 0; i < m_cube.size(); ++i) {
            m_cube[i]->counts(0.0);
        }

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill model into model cube
 *
 * @param[in] obs CTA observation.
 *
 * Adds the expected number of events for a given observation to the events
 * that are already found in the model cube. The method also updates the
 * GTI of the model cube so that cube GTI is a list of the GTIs of all
 * observations that were used to generate the model cube.
 ***************************************************************************/
void ctmodel::fill_cube(const GCTAObservation* obs)
{
    // Continue only if observation pointer is valid
    if (obs != NULL) {

        // Get energy boundaries and GTI references for observation
        const GEbounds& ebounds = obs->events()->ebounds();
        const GGti&     gti     = obs->events()->gti();

        // Initialise statistics
        double sum              = 0.0;
        int    num_outside_ebds = 0;

        // Setup cube GTIs for this observation
        m_cube.gti(obs->events()->gti());

        // Loop over all cube bins
        for (int i = 0; i < m_cube.size(); ++i) {

            // Get cube bin
            GCTAEventBin* bin = m_cube[i];

            // Skip bin if it is outside the energy range of the observation
            if (!ebounds.contains(bin->energy())) {
                num_outside_ebds++;
                continue;
            }

            // Get actual bin value
            double value = bin->counts();
            
            // Compute model value for cube bin
            double model = m_obs.models().eval(*bin, *obs) * bin->size();

            // Add model to actual value
            value += model;
            sum   += model;

            // Store value
            bin->counts(value);

        } // endfor: looped over all cube bins

        // Append GTIs of observation to list of GTIs
        m_gti.extend(gti);

        // Update GTIs
        m_cube.gti(m_gti);

        // Log results
        if (logTerse()) {
            log << gammalib::parformat("Model events in cube");
            log << sum << std::endl;
            log << gammalib::parformat("Cube bins outside energy range");
            log << num_outside_ebds << std::endl;
        }

        // Log cube
        if (logExplicit()) {
            log.header2("Model cube");
            log << m_cube << std::endl;
        }

    } // endif: observation was valid

    // Return
    return;
}
