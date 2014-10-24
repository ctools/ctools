/***************************************************************************
 *                        ctool - ctool base class                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @brief Get the energy boundaries
 *
 * @exception GException::invalid_value
 *            No valid energy boundary extension found.
 *
 * Get the energy boundaries according to the user parameters. The method
 * supports loading of energy boundary information from the EBOUNDS or
 * ENERGYBINS extension, or setting energy boundaries using a linear
 * or logarithmical spacing.
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
 * @brief Set response for all CTA observations in container
 *
 * @param[in,out] obs Observation container
 *
 * Set the response for a CTA observations in the container that so far have
 * no response using the "database" and "irf" task parameters.
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
 * Set the response for a CTA observation using the "database" and "irf"
 * task parameters.
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
