/***************************************************************************
 *             ctobservation - Base class for observation tools            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2017 by Juergen Knoedlseder                         *
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
 * @file ctobservation.cpp
 * @brief Observation tool base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ctobservation.hpp"

/* __ Method name definitions ____________________________________________ */

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
 * @param[in] name Observation tool name.
 * @param[in] version Observation tool version.
 *
 * Constructs a observation tool from the @p name and @p version. See the
 * equivalent ctool constructor for details.
 ***************************************************************************/
ctobservation::ctobservation(const std::string& name,
                             const std::string& version) : ctool(name, version)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observations constructor
 *
 * @param[in] name Observation tool name.
 * @param[in] version Observation tool version.
 * param[in] obs Observation container.
 *
 * Constructs a observation tool from the @p name, @p version and an
 * observation container.
 ***************************************************************************/
ctobservation::ctobservation(const std::string&   name,
                             const std::string&   version,
                             const GObservations& obs) : ctool(name, version)
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
 * @param[in] name Observation tool name.
 * @param[in] version Observation tool version.
 * @param[in] argc Number of arguments in command line.
 * @param[in] argv Array of command line arguments.
 *
 * Constructs a observation tool from the @p name, @p version and command
 * line arguments. See the equivalent ctool constructor for details.
 ***************************************************************************/
ctobservation::ctobservation(const std::string& name,
                             const std::string& version,
                             int   argc,
                             char *argv[]) : ctool(name, version, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app Observation tool.
 *
 * Constructs an instance of a observation tool by copying information from
 * another observation tool.
 ***************************************************************************/
ctobservation::ctobservation(const ctobservation& app) : ctool(app)
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
 *
 * Destructs the observation tool.
 ***************************************************************************/
ctobservation::~ctobservation(void)
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
 * @param[in] app Observation tool.
 * @return Observation tool.
 *
 * Assigns a observation tool.
 ***************************************************************************/
ctobservation& ctobservation::operator=(const ctobservation& app)
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
 * @brief Write OGIP keywords in FITS HDU
 *
 * @param[in,out] hdu Pointer to FITS HDU.
 *
 * Writes OGIP keywords in FITS HDU.
 ***************************************************************************/
void ctobservation::write_ogip_keywords(GFitsHDU* hdu) const
{
    // Continue only if pointer is valid
    if (hdu != NULL) {

        // Set creator
        std::string creator = this->name() + " v" + this->version();

        // Compute times
        std::string utc_obs  = m_ogip_tstart.utc();
        std::string utc_end  = m_ogip_tstop.utc();
        std::string date_obs = utc_obs.substr(0, 10);
        std::string time_obs = utc_obs.substr(11, 8);
        std::string date_end = utc_end.substr(0, 10);
        std::string time_end = utc_end.substr(11, 8);

        // Compute deadtime correction
        double deadc = (m_ogip_ontime > 0.0) ? m_ogip_livetime / m_ogip_ontime : 1.0;

        // Set observation information
        hdu->card("CREATOR",  creator,          "Program which created the file");
        hdu->card("TELESCOP", m_ogip_telescope, "Telescope");

        // Set observation time information
        hdu->card("DATE_OBS", date_obs,        "Observation start date");
        hdu->card("TIME_OBS", time_obs,        "Observation start time");
        hdu->card("DATE_END", date_end,        "Observation end date");
        hdu->card("TIME_END", time_end,        "Observation end time");
        hdu->card("TELAPSE",  m_ogip_telapse,  "[s] Elapsed time");
        hdu->card("ONTIME",   m_ogip_ontime,   "[s] Total good time including deadtime");
        hdu->card("LIVETIME", m_ogip_livetime, "[s] Total livetime");
        hdu->card("EXPOSURE", m_ogip_exposure, "[s] Exposure time");
        hdu->card("DEADC",    deadc,           "Deadtime correction factor");

    } // endif: pointer was valid

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
void ctobservation::init_members(void)
{
    // Initialise protected members
    m_obs.clear();

    // Initialise private members
    m_ogip_telescope.clear();
    m_ogip_tstart.clear();
    m_ogip_tstop.clear();
    m_ogip_telapse   = 0.0;
    m_ogip_exposure  = 0.0;
    m_ogip_ontime    = 0.0;
    m_ogip_livetime  = 0.0;
    m_index_unbinned = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app Observation tool.
 ***************************************************************************/
void ctobservation::copy_members(const ctobservation& app)
{
    // Copy protected members
    m_obs = app.m_obs;

    // Copy private members
    m_ogip_telescope = app.m_ogip_telescope;
    m_ogip_tstart    = app.m_ogip_tstart;
    m_ogip_tstop     = app.m_ogip_tstop;
    m_ogip_telapse   = app.m_ogip_telapse;
    m_ogip_exposure  = app.m_ogip_exposure;
    m_ogip_ontime    = app.m_ogip_ontime;
    m_ogip_livetime  = app.m_ogip_livetime;
    m_index_unbinned = app.m_index_unbinned;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctobservation::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return first unbinned CTA observation (const version)
 *
 * @return Const pointer to first unbinned CTA observation
 *
 * Returns a const pointer to the first unbinned CTA observation in the
 * container. If no CTA observation exists a NULL pointer is returned.
 *
 * The method calls next_unbinned_observation(). See the method for details.
 ***************************************************************************/
const GCTAObservation* ctobservation::first_unbinned_observation(void) const
{
    // Initialise index
    m_index_unbinned = 0;

    // Initialise OGIP members
    m_ogip_telescope.clear();
    m_ogip_tstart.clear();
    m_ogip_tstop.clear();
    m_ogip_telapse  = 0.0;
    m_ogip_exposure = 0.0;
    m_ogip_ontime   = 0.0;
    m_ogip_livetime = 0.0;

    // Get next unbinned CTA observation
    const GCTAObservation* obs = next_unbinned_observation();

    // Return first CTA observation
    return obs;
}


/***********************************************************************//**
 * @brief Return next unbinned CTA observation (const version)
 *
 * @return Const pointer to next unbinned CTA observation
 *
 * Returns a const pointer to the next unbinned CTA observation in the
 * container. If no CTA observation exists any more a NULL pointer is
 * returned.
 *
 * The method writes for each encountered observation a level 3 header into
 * the logger. It will also signal when an observation was skipped because
 * it either was not a CTA observation or not an unbinned observation.
 *
 * @todo Logger methods should be declared const to avoid the const casting.
 ***************************************************************************/
const GCTAObservation* ctobservation::next_unbinned_observation(void) const
{
    // Initialise pointer on CTA observation
    const GCTAObservation* obs = NULL;

    // Loop over all remaining observation in the container
    for (; m_index_unbinned < m_obs.size(); ++m_index_unbinned) {

        // Write header for the current observation
        const_cast<ctobservation*>(this)->log_header3(TERSE,
                                   get_obs_header(m_obs[m_index_unbinned]));

        // Case the observation to a CTA observation. This will return a
        // NULL pointer if the observation is not a CTA observation.
        obs = dynamic_cast<const GCTAObservation*>(m_obs[m_index_unbinned]);

        // Skip observation if it's not CTA
        if (obs == NULL) {
            std::string msg = " Skipping "+
                              m_obs[m_index_unbinned]->instrument()+
                              " observation";
            const_cast<ctobservation*>(this)->log_string(NORMAL, msg);
            continue;
        }

        // Skip observation if we have a binned observation
        if (obs->eventtype() == "CountsCube") {
            obs             = NULL;
            std::string msg = " Skipping binned "+
                              m_obs[m_index_unbinned]->instrument()+
                              " observation";
            const_cast<ctobservation*>(this)->log_string(NORMAL, msg);
            continue;
        }

        // If we come to this point we have an unbinned CTA observation and
        // we can forward the index to the next index and break the loop
        m_index_unbinned++;
        break;

    } // endfor: looped over all observation

    // Update OGIP members
    if (obs != NULL) {
        if (m_ogip_telapse == 0.0) {
            m_ogip_telescope = obs->instrument();
            m_ogip_tstart    = obs->gti().tstart();
            m_ogip_tstop     = obs->gti().tstop();
            m_ogip_telapse   = obs->gti().telapse();
            m_ogip_exposure  = obs->ontime();
            m_ogip_ontime    = obs->ontime();
            m_ogip_livetime  = obs->livetime();
        }
        else {
            if (obs->gti().tstart() < m_ogip_tstart) {
                m_ogip_tstart = obs->gti().tstart();
            }
            if (obs->gti().tstop() > m_ogip_tstop) {
                m_ogip_tstop = obs->gti().tstop();
            }
            m_ogip_telapse  += obs->gti().telapse();
            m_ogip_exposure += obs->ontime();
            m_ogip_ontime   += obs->ontime();
            m_ogip_livetime += obs->livetime();
        }
    }

    // Return next CTA observation
    return obs;
}


/***********************************************************************//**
 * @brief Save event list in FITS format.
 *
 * Save the event list as a FITS file. The file name of the FITS file is
 * specified by the `outobs` parameter.
 ***************************************************************************/
void ctobservation::save_events_fits(void)
{
    // Save only if there are observations
    if (m_obs.size() > 0) {

        // Get output filename
        std::string outobs = (*this)["outobs"].filename();

        // Loop over all unbinned CTA observations in the container
        for (GCTAObservation* obs = first_unbinned_observation(); obs != NULL;
             obs = next_unbinned_observation()) {

            // Get input file name and default extension names
            std::string infile  = obs->eventfile();
            std::string evtname = gammalib::extname_cta_events;
            std::string gtiname = gammalib::extname_gti;

            // Create file name object
            GFilename fname(outobs);

            // Extract filename and event extension name
            std::string outfile = fname.url();

            // Append event extension name. We handle here the possibility
            // to write the events into a different extension.
            if (fname.has_extname()) {
                outfile += "["+fname.extname()+"]";
            }
            else {
                outfile += "["+gammalib::extname_cta_events+"]";
            }

            // Log filename
            log_value(NORMAL, "Event list file", outfile);
            
            // Save event list
            save_event_list(obs, infile, evtname, gtiname, outfile);

            // Exit the loop (if this method is called we should only
            // have a single unbinned observation in the observation
            // container
            break;

        } // endfor: looped over all unbinned observations

    } // endif: there were observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save event list(s) in XML format.
 *
 * Save the event list(s) into FITS files and write the file path information
 * into a XML file. The filename of the XML file is specified by the `outobs`
 * parameter, the filename(s) of the event lists are built by prepending a
 * prefix to the input event list filenames. Any path present in the input
 * filename will be stripped, i.e. the event list(s) will be written in the
 * local working directory (unless a path is specified in the prefix).
 ***************************************************************************/
void ctobservation::save_events_xml(void)
{
    // Get output filename and prefix
    std::string outobs = (*this)["outobs"].filename();

    // Issue warning if output filename has no .xml suffix
    log_string(TERSE, warn_xml_suffix(outobs));

    // Loop over all unbinned CTA observations in the container
    for (GCTAObservation* obs = first_unbinned_observation(); obs != NULL;
         obs = next_unbinned_observation()) {

        // Get input file name and default extension names
        std::string infile  = obs->eventfile();
        std::string evtname = gammalib::extname_cta_events;
        std::string gtiname = gammalib::extname_gti;

        // Extract event and GTI extension names from input FITS file
        GFilename fname(infile);
        if (fname.has_extname()) {
            evtname = fname.extname();
        }
        gtiname = get_gtiname(fname.url(), evtname);

        // Set event output file name
        std::string outfile = set_outfile_name(infile);

        // Append event extension name
        outfile += "["+evtname+"]";

        // Log filename
        log_value(NORMAL, "Event list file", outfile);

        // Store output file name in observation
        obs->eventfile(outfile);

        // Save event list
        save_event_list(obs, infile, evtname, gtiname, outfile);

    } // endfor: looped over observations

    // Write observation definition XML file name into logger
    log_value(NORMAL, "Obs. definition file", outobs);

    // Save observations in XML file
    m_obs.save(outobs);

    // Return
    return;
}
