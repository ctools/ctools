/***************************************************************************
 *                  ctmapcube - Map cube generation tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file ctmapcube.cpp
 * @brief Map cube generation tool implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "ctmapcube.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
ctmapcube::ctmapcube(void) : ctool(CTMAPCUBE_NAME, CTMAPCUBE_VERSION)
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
ctmapcube::ctmapcube(int argc, char *argv[]) :
           ctool(CTMAPCUBE_NAME, CTMAPCUBE_VERSION, argc, argv)
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
ctmapcube::ctmapcube(const ctmapcube& app) : ctool(app)
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
ctmapcube::~ctmapcube(void)
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
ctmapcube& ctmapcube::operator=(const ctmapcube& app)
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
void ctmapcube::clear(void)
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
 * @brief Generate map cube.
 *
 * This method reads the task parameters from the parfile and generates the
 * map cube.
 ***************************************************************************/
void ctmapcube::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }

    // Get task parameters
    get_parameters();

    // Write models into logger
    if (logTerse()) {
        log << std::endl;
        log.header1(gammalib::number("Model", m_models.size()));
        log << m_models << std::endl;
    }

    // Generate map cube
    create_cube();

    // Optionally publish map cube
    if (m_publish) {
        publish();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save map cube
 *
 * Saves the map cube into a FITS file.
 ***************************************************************************/
void ctmapcube::save(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Save map cube");
    }

    // Get map cube filename
    m_outcube = (*this)["outcube"].filename();

    // Save only if filename is non-empty
    if (!m_outcube.is_empty()) {

        // Log filename
        if (logTerse()) {
            log << gammalib::parformat("Map cube file");
            log << m_outcube.url() << std::endl;
        }

        // Save map cube into FITS file
        m_cube.save(m_outcube, clobber());

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Publish map cube
 *
 * @param[in] name Map cube name.
 ***************************************************************************/
void ctmapcube::publish(const std::string& name)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Publish map cube");
    }

    // Set default name is user name is empty
    std::string user_name(name);
    if (user_name.empty()) {
        user_name = CTMAPCUBE_NAME;
    }

    // Log filename
    if (logTerse()) {
        log << gammalib::parformat("Map cube");
        log << user_name << std::endl;
    }

    // Publish map cube
    m_cube.cube().publish(user_name);

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
void ctmapcube::init_members(void)
{
    // Initialise members
    m_outcube.clear();
    m_publish = false;

    // Initialise protected members
    m_models.clear();
    m_cube.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Application.
 ***************************************************************************/
void ctmapcube::copy_members(const ctmapcube& app)
{
    // Copy attributes
    m_outcube = app.m_outcube;
    m_publish = app.m_publish;

    // Copy protected members
    m_models = app.m_models;
    m_cube   = app.m_cube;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctmapcube::free_members(void)
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
void ctmapcube::get_parameters(void)
{
    // Get energy binning parameters
    GEnergies energies = create_energies();

    // Get spatial binning parameters
    double      xref     = (*this)["xref"].real();
    double      yref     = (*this)["yref"].real();
    std::string proj     = (*this)["proj"].string();
    std::string coordsys = (*this)["coordsys"].string();
    double      binsz    = (*this)["binsz"].real();
    int         nxpix    = (*this)["nxpix"].integer();
    int         nypix    = (*this)["nypix"].integer();

    // Read model definition file if required
    if (m_models.size() == 0) {

        // Get model filename
        GFilename inmodel = (*this)["inmodel"].filename();

        // Load models from file
        m_models.load(inmodel);

    } // endif: there were no models

    // Get remaining parameters
    m_publish = (*this)["publish"].boolean();

    // Read optionally output cube filenames
    if (read_ahead()) {
        m_outcube = (*this)["outcube"].filename();
    }

    // Allocate map cube
    GSkyMap cube(proj, coordsys, xref, yref, -binsz, binsz, nxpix, nypix,
                 energies.size());

    // Allocate map cube
    m_cube = GModelSpatialDiffuseCube(cube, energies);

    // Write parameters into logger
    if (logTerse()) {
        log_parameters();
        log << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create energy vector from user parameters
 *
 * Get the energy vector according to the user parameters. The method
 * supports loading of energy information from the ENERGIES extension (or
 * any other extension specified by the "ebinfile" parameters), or setting
 * energies using a linear or logarithmical spacing.
 *
 * The following parameters are read:
 *      ebinalg - Energy binning algorithm
 *      emin - Minimum energy (if ebinalg != FILE)
 *      emax - Maximum energy (if ebinalg != FILE)
 *      enumbins - Number of energies (if ebinalg != FILE)
 ***************************************************************************/
GEnergies ctmapcube::create_energies(void)
{
    // Allocate energies
    GEnergies energies;

    // Get energy binning algorithm
    std::string ebinalg = (*this)["ebinalg"].string();

    // If energy binning algorithm is of type "FILE" (case sensitive), then
    // read energies from FITS file ...
    if (ebinalg == "FILE") {

        // Get filename
        GFilename ebinfile = (*this)["ebinfile"].filename();

        // Load energies
        energies.load(ebinfile);

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

        // Setup energies
        energies = GEnergies(enumbins, GEnergy(emin, "TeV"),
                                       GEnergy(emax, "TeV"), log);

    } // endelse: ebinalg was not "FILE"

    // Return energies
    return energies;
}


/***********************************************************************//**
 * @brief Generate map cube
 *
 * Generate map cube by looping over all sky models in the container.
 ***************************************************************************/
void ctmapcube::create_cube(void)
{
    // Write header
    if (logTerse()) {
        log << std::endl;
        log.header1("Generate map cube");
    }

    // Loop over all models
    for (int i = 0; i < m_models.size(); ++i) {

        // Get pointer to sky model
        GModelSky* model = dynamic_cast<GModelSky*>(m_models[i]);

        // Fall through if model is not a sky model
        if (model == NULL) {

            // Signal that model is skipped
            if (logNormal()) {
                log << gammalib::parformat("Skip model");
                log << m_models[i]->name() << std::endl;
            }

            // Continue
            continue;

        }

        // Signal that model is used
        if (logNormal()) {
            log << gammalib::parformat("Use model");
            log << m_models[i]->name() << std::endl;
        }

        // Add sky model
        if (dynamic_cast<GModelSpatialPointSource*>(model->spatial()) != NULL) {
            add_ptsrc_model(model);
        }
        else {
            add_model(model);
        }

    } // endfor: looped over energies

    // Log cube
    if (logTerse()) {
        log << std::endl;
        log.header1("Map cube");
        log << m_cube << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Add one model to map cube
 *
 * Adds one model to map cube.
 ***************************************************************************/
void ctmapcube::add_model(GModelSky* model)
{
    // Get cube energies
    GEnergies energies = m_cube.energies();

    // Get pointer to sky map
    GSkyMap* map = const_cast<GSkyMap*>(&m_cube.cube());

    // Loop over all energy bins
    for (int iebin = 0; iebin < energies.size(); ++iebin) {

        // Loop over all spatial pixels
        for (int ipixel = 0; ipixel < map->npix(); ++ipixel) {

            // Set photon
            GPhoton photon(map->inx2dir(ipixel), energies[iebin], GTime());

            // Add model value
            map->operator()(ipixel, iebin) += model->value(photon);

        } // endfor: looped over all pixels

    } // endfor: looped over all energies

    // Return
    return;
}


/***********************************************************************//**
 * @brief Add one point source model to map cube
 *
 * Adds one point source model to map cube.
 ***************************************************************************/
void ctmapcube::add_ptsrc_model(GModelSky* model)
{
    // Get point source spatial component
    GModelSpatialPointSource* ptsrc = static_cast<GModelSpatialPointSource*>(model->spatial());

    // Get point source direction
    GSkyDir dir = ptsrc->dir();

    // Get cube energies
    GEnergies energies = m_cube.energies();

    // Get pointer to sky map
    GSkyMap* map = const_cast<GSkyMap*>(&m_cube.cube());

    // Get test sky map for flux computation
    GSkyMap test_map(*map);

    // Continue only if map contains point source sky direction
    if (map->contains(dir)) {

        // Get map pixel
        int ipixel = map->dir2inx(dir);

        // Loop over all energy bins
        for (int iebin = 0; iebin < energies.size(); ++iebin) {

            // Set photon using the point source sky direction
            GPhoton photon(dir, energies[iebin], GTime());

            // Compute model value
            double value = model->value(photon);

            // Compute flux normalization
            test_map(ipixel, iebin) = value;
            double flux = test_map.flux(ipixel, iebin);
            double norm = (flux > 0.0) ? value / flux : 0.0;

            // Add model to sky map
            map->operator()(ipixel, iebin) += norm * value;

        } // endfor: looped over all energies

    } // endif: map contained sky direction

    // Return
    return;
}
