/***************************************************************************
 *                    ctfilter - CTA cube filter tool                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Chia-Chun Lu                                     *
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
 * @brief CTA mask cube implementation
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
#define G_RUN                                               "ctcubemask::run()"
#define G_APPLY_MASK                "ctcubemask::apply_mask(GCTAObservation*,"\
                                                              "std::string&)"

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
ctcubemask::ctcubemask(void) : GApplication(CTCUBEMASK_NAME, CTCUBEMASK_VERSION)
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
 * This constructor creates an instance of the class that is initialised from
 * an observation container.
 ***************************************************************************/
ctcubemask::ctcubemask(GObservations obs) : GApplication(CTCUBEMASK_NAME, CTCUBEMASK_VERSION)
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
ctcubemask::ctcubemask(int argc, char *argv[]) : 
                   GApplication(CTCUBEMASK_NAME, CTCUBEMASK_VERSION, argc, argv)
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
ctcubemask::ctcubemask(const ctcubemask& app) : GApplication(app)
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
 ***************************************************************************/
ctcubemask& ctcubemask::operator= (const ctcubemask& app)
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
void ctcubemask::clear(void)
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
 * This method performs the event selection and saves the result
 ***************************************************************************/
void ctcubemask::execute(void)
{
    // Signal that some parameters should be read ahead
    m_read_ahead = true;
    
    // Perform event selection
    run();

    // Save results
    save();

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
        log.header1("Observations for applying mask");
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

            // Get temporary file name
            std::string filename = std::tmpnam(NULL);

            // Save observation in temporary file
            obs->save(filename, true);

            // Log saved FITS file
            if (logExplicit()) {
                GFits tmpfile(filename);
                log << std::endl;
                log.header1("FITS file content of temporary file");
                log << tmpfile << std::endl;
                tmpfile.close();
            }

            // Check temporary file
            std::string message = check_infile(filename);
            if (message.length() > 0) {
                throw GException::app_error(G_RUN, message);
            }

            // Load observation from temporary file, including event selection
            apply_mask(obs, filename);

            // Remove temporary file
            std::remove(filename.c_str());
            
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
        log.header1("Observations after selection");
        log << m_obs << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save the selected event list(s)
 *
 * This method saves the selected event list(s) into FITS file(s). There are
 * two modes, depending on the m_use_xml flag.
 *
 * If m_use_xml is true, all selected event list(s) will be saved into FITS
 * files, where the output filenames are constructued from the input
 * filenames by prepending the m_prefix string to name. Any path information
 * will be stripped form the input name, hence event files will be written
 * into the local working directory (unless some path information is present
 * in the prefix). In addition, an XML file will be created that gathers
 * the filename information for the selected event list(s). If an XML file
 * was present on input, all metadata information will be copied from this
 * input file.
 *
 * If m_use_xml is false, the selected event list will be saved into a FITS
 * file.
 ***************************************************************************/
void ctcubemask::save(void)
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

    // Case A: Save event file(s) and XML metadata information
    if (m_use_xml) {
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
 * @brief Get application parameters
 *
 * Get all task parameters from parameter file or (if required) by querying
 * the user. Times are assumed to be in the native CTA MJD format.
 *
 * This method also loads observations if no observations are yet allocated.
 * Observations are either loaded from a single CTA even list, or from a
 * XML file using the metadata information that is stored in that file.
 ***************************************************************************/
void ctcubemask::get_parameters(void)
{
    // If there are no observations in container then add a single CTA
    // observation using the parameters from the parameter file
    if (m_obs.size() == 0) {

        // Get CTA event list file name
        m_infile = (*this)["infile"].filename();

        // Allocate CTA observation
        GCTAObservation obs;

        // Try first to open as FITS file
        try {

            // Load event list in CTA observation
            obs.load_binned(m_infile);

            // Append CTA observation to container
            m_obs.append(obs);

            // Signal that no XML file should be used for storage
            m_use_xml = false;
            
        }
        
        // ... otherwise try to open as XML file
        catch (GException::fits_open_error &e) {

            // Load observations from XML file
            m_obs.load(m_infile);

            // Signal that XML file should be used for storage
            m_use_xml = true;

        }

    } // endif: there was no observation in the container

    // Get parameters
	m_regfile = (*this)["regfile"].filename();
    m_usepnt = (*this)["usepnt"].boolean();
    if (!m_usepnt) {
        m_ra  = (*this)["ra"].real();
        m_dec = (*this)["dec"].real();
    }
    m_rad  = (*this)["rad"].real();
    m_emin = (*this)["emin"].real();
    m_emax = (*this)["emax"].real();
	
    // Optionally read ahead parameters so that they get correctly
    // dumped into the log file
    if (m_read_ahead) {
        m_outfile = (*this)["outfile"].filename();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Select events
 *
 * @param[in] obs CTA observation.
 * @param[in] filename File name.
 *
 * Select events from a FITS file by making use of the selection possibility
 * of the cfitsio library on loading a file. A selection string is created
 * from the specified criteria that is appended to the filename so that
 * cfitsio will automatically filter the event data. This selection string
 * is then applied when opening the FITS file. The opened FITS file is then
 * saved into a temporary file which is the loaded into the actual CTA
 * observation, overwriting the old CTA observation. The ROI, GTI and EBounds
 * of the CTA event list are then set accordingly to the specified selection.
 * Finally, the temporary file created during this process is removed.
 *
 * Good Time Intervals of the observation will be limited to the time
 * interval [m_tmin, m_tmax]. If m_tmin=m_tmax=0, no time selection is
 * performed.
 *
 * @todo Use INDEF instead of 0.0 for pointing as RA/DEC selection
 ***************************************************************************/
void ctcubemask::apply_mask(GCTAObservation* obs, const std::string& filename)
{
  std::cout<< "start appply filter" << std::endl;
    // Open FITS file
    GFits file(filename);
    GCTAEventCube* cube = (GCTAEventCube*) obs->events();
    GEbounds ebounds = cube->ebounds();
    GSkymap map = cube->map();

    // Apply energy selection
    int e_idx1 = 0;
    int e_idx2 = 0;
    int n_ebin = ebounds.size();

    // Loop over ebounds to find the first energy band
    for ( int i = 0 ; i < ebounds.size() ; i++){
      double emin = ebounds.emin(i).TeV();
      if ( emin > m_emin ){
	e_idx1 = i;
	break;
      }
    }

    std::cout << "finish looping over ebounds to find 1st" << std::endl;
    std::cout << "e_ix1 " << e_idx1 << std::endl;

    // Loop over ebounds to find the last energy band
    for ( int i = ebounds.size()-1 ; i < ebounds.size() ; i--){
      double emin = ebounds.emin(i).TeV();
      if ( emin < m_emax ){
	e_idx2 = i-1;
	break;
      }
    }
     std::cout << "finish looping over ebounds to find 1st" << std::endl;
     std::cout << "e_ix2 " << e_idx2 << std::endl;

    // Set all pixels outside the desired energy bands negative
    int npix = map.npix();
    for ( int i = 0 ; i < e_idx1 ; i++){
      for ( int pixel = 0 ; pixel < npix ; pixel++){
		  map(pixel,i) = -1;
      }
    }
    for ( int i = e_idx2 ; i < ebounds.size() ; i++){
      for ( int pixel = 0 ; pixel < npix ; pixel++){
		  map(pixel,i) = -1;
      }
    }
    std::cout << "finish energy band filter" << std::endl;


   
    // Set all pixels inside the desired energy bands 
    // but outside ROI or inside exlusion regions negative
    GSkyRegions regs = GSkyRegions(m_regfile);
	GSkyRegionCircle roi = GSkyRegionCircle(m_ra, m_dec, m_rad);
    for ( int i = e_idx1 ; i <= e_idx2 ; i++){
      for ( int pixel = 0 ; pixel < npix ; pixel++){
		  GSkyDir dir = map.inx2dir(pixel);
		  if ( roi.contains(dir) == false || 
			   regs.contains(dir) == true){
			  map(pixel,i) = -1;
		  }
      }
    }
    std::cout << "finish roi selection" << std::endl;
   
    GCTAEventCube newcube( map, ebounds, obs->events()->gti());
    obs->events(newcube);

    // Log selected FITS file
    if (logExplicit()) {
        log << std::endl;
        log.header1("FITS file content after filter");
        log << file << std::endl;
    }

    // Check if we have an IMAGE HDU
    if (!file.contains("IMAGE")) {
        std::string message = "No \"IMAGE\" extension found in FITS file. ";
        throw GException::app_error(G_APPLY_MASK, message);
    }


    // Get CTA event list pointer
    GCTAEventList* list =
        static_cast<GCTAEventList*>(const_cast<GEvents*>(obs->events()));

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
    m_infile.clear();
	m_regfile.clear();
    m_outfile.clear();
	m_prefix.clear();
    m_usepnt = false;
    m_ra     = 0.0;
    m_dec    = 0.0;
    m_rad    = 0.0;
    m_emin   = 0.0;
    m_emax   = 0.0;

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
void ctcubemask::copy_members(const ctcubemask& app)
{
    // Copy parameters
    m_infile  = app.m_infile;
	m_regfile = app.m_regfile;
    m_outfile = app.m_outfile;
	m_prefix  = app.m_prefix;
    m_usepnt  = app.m_usepnt;
    m_ra      = app.m_ra;
    m_dec     = app.m_dec;
    m_rad     = app.m_rad;
    m_emin    = app.m_emin;
    m_emax    = app.m_emax;

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
void ctcubemask::free_members(void)
{
    // Write separator into logger
    if (logTerse()) {
        log << std::endl;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check input filename
 *
 * @param[in] filename File name.
 *
 * This method checks if the input FITS file is correct.
 ***************************************************************************/
std::string ctcubemask::check_infile(const std::string& filename) const
{
    // Initialise message string
    std::string message = "";

    // Open FITS file
    GFits file(filename);

    // Check for IMAGE HDU
    GFitsImage* image = NULL;
    try {

        // Get pointer to FITS image
        image = file.image("IMAGE");
    }
    catch (GException::fits_hdu_not_found& e) {
        message = "No \"IMAGE\" extension found in input file \""
                + m_outfile + "\".";
    }

    // Return
    return message;
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
 * @brief Save counts map in FITS format.
 *
 * Save the counts map as a FITS file. The filename of the FITS file is
 * specified by the m_outfile member.
 ***************************************************************************/
void ctcubemask::save_fits(void)
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
void ctcubemask::save_xml(void)
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

