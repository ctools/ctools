/***************************************************************************
 *   ctfindvar - search time variability tool                              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Simon Bonnefoy                                   *
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
 * @file ctfindvar.cpp
 * @brief search time variability tool implementation
 * @author Simon Bonnefoy
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ctfindvar.hpp"
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include <cmath>
#include <sstream>

/* __ OpenMP section _____________________________________________________ */
//#ifdef _OPENMP
//#include <omp.h>
//#endif

/* __ Method name definitions ____________________________________________ */
#define G_FILL_CUBE                  "ctfindvar::fill_cube(GCTAObservation*)"
#define G_INIT_CUBE                              "ctfindvar::init_cube(void)"
#define G_INIT_GTIS                              "ctfindvar::init_gtis(void)"

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
 *
 * Constructs empty search time variability tool.
 ***************************************************************************/
ctfindvar::ctfindvar(void) :
    ctobservation(CTFINDVAR_NAME, VERSION)
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
 * Constructs search time variability tool from an observation container.
 ***************************************************************************/
ctfindvar::ctfindvar(const GObservations& obs) :
    ctobservation(CTFINDVAR_NAME, VERSION, obs)
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
 *
 * Constructs search time variability tool using command line arguments for user
 * parameter setting.
 ***************************************************************************/
ctfindvar::ctfindvar(int argc, char *argv[]) :
    ctobservation(CTFINDVAR_NAME, VERSION, argc, argv)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] app search time variability tool.
 *
 * Constructs search time variability tool from another search time variability tool.
 ***************************************************************************/
ctfindvar::ctfindvar(const ctfindvar& app) :
    ctobservation(app)
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
 * Destructs search time variability tool.
 ***************************************************************************/
ctfindvar::~ctfindvar(void)
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
 * @param[in] app search time variability tool.
 * @return search time variability tool.
 *
 * Assigns search time variability tool.
 ***************************************************************************/
ctfindvar& ctfindvar::operator=(const ctfindvar& app)
{
    // Execute only if object is not identical
    if (this != &app) {

        // Copy base class members
        this->ctobservation::operator=(app);

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
 * @brief Clear search time variability tool
 *
 * Clears search time variability tool.
 ***************************************************************************/
void ctfindvar::clear(void)
{
    // Free members
    free_members();
    this->ctobservation::free_members();
    this->ctool::free_members();

    // Clear base class (needed to conserve tool name and version)
    this->GApplication::clear();

    // Initialise members
    this->ctool::init_members();
    this->ctobservation::init_members();
    init_members();

    // Write header into logger
    log_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Run search time variability tool
 ***************************************************************************/
void ctfindvar::run(void)
{
    // If we're in debug mode then all output is also dumped on the screen
    if (logDebug()) {
        log.cout(true);
    }
    
    // Get task parameters
    get_parameters();

    // Write input observation container into logger
    log_observations(NORMAL, m_obs, "Input observation");

     // Loop over all unbinned CTA observations in the container
    #pragma omp parallel for
    for (int i = 0; i < m_obs.size(); ++i) {

        // Get pointer to observation
        GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[i]);

        // Fill the cube
        fill_cube(obs);

        // Dispose events to free memory
        obs->dispose_events();

    } // endfor: looped over observations

    // Smooth the maps if requested
    if ((*this)["smoothkrnl"].is_valid() && 
        (*this)["smoothkrnl"].string() != "NONE") {
        m_counts.smooth((*this)["smoothkrnl"].string(),
                        (*this)["smoothpar"].real());
    }

    const int nbins = m_counts.nmaps();

    //preparing the histogram with significance evolution for the source center and the highest sig-pix
    GNdarray pixSig(nbins), pixSigSrc(nbins), pixSigMax(nbins); //Storing the significance for each GTI of the pixels
  
    // creating GSkyDir to get the position of the source and each pixels
    GSkyDir srcSkyDir;

    if ((*this)["coordsys"].string()=="CEL")
    {
        srcSkyDir.radec_deg((*this)["xsrc"].real(),(*this)["ysrc"].real());
    }
    else
    {
        srcSkyDir.lb_deg((*this)["xsrc"].real(),(*this)["ysrc"].real());
    }
    //Getting the index of pixel on which falls the source of interest
    int srcInxPix = m_counts.dir2inx(srcSkyDir);

    //Prepare the final skymap with the max significance of each pixel
    GSkyMap pixVarSig = m_counts.extract(1);
    double max_sig=0;

    //looping over all the pixels in the cube
    for (int pix_number=0; pix_number<m_counts.npix(); pix_number++)
    {    
        double total_counts=0;
        for (int i=0;i<m_counts.nmaps();i++)
        {
            total_counts += m_counts(pix_number, i);
        }
        if(pix_number%20==0)
        {
            std::cout << "checking pixel number " << pix_number << std::endl;
            std::cout << "number of counts in pixel" << total_counts<< std::endl;
        }

        get_variability_sig(pix_number,nbins, pixSig);
        if (srcInxPix == pix_number) //Getting the significance evolution for the source
        {
            pixSigSrc=pixSig;
            std::cout << "checking pixel number of the source of interest" << pix_number << std::endl;
            std::cout << "number of counts in pixel of the source of interest" << total_counts << std::endl;
            std::cin.ignore();
        }
        else
        {
            // Getting the evolution for the pix with highest significance
            if (max(pixSig) > max_sig)  
            {
                pixSigMax=pixSig;
            }
        }
        //storing sig in skymap
        pixVarSig(pix_number) = max(pixSig);
    }
    return;
}


/***********************************************************************//**
 * @brief Obtain the significance of variability in a given bin
 * 
 * @param[in]  pix_number        Pixel in the skymap
 * @param[in]  nbins             Number of bins in the skymap
 * @param[out] sig_histogram     Histogram for storing returned significances
 ***************************************************************************/
void ctfindvar::get_variability_sig(const int& pix_number, 
                                    const int& nbins, 
                                    GNdarray&  sig_histogram)
{
    std::vector<bool> accepted_bin_bckg_vector;
    std::vector<double> excess_bin_vector;
    int background_bin_array[nbins];
    bool background_validated=false;
    int non, noff;
    double alpha, sig;

    for (int i=0;i<nbins;i++)
    {
        background_bin_array[i]=0;
        accepted_bin_bckg_vector.push_back(1);
        excess_bin_vector.push_back(0);
    }

    while (background_validated==false)
    {
        background_validated=true;

        for (int i=0; i< nbins; i++) //looping over all the GTIs of the pixel
        {
           alpha=0;
           if (accepted_bin_bckg_vector[i]==0) continue;     //the GTI is discared from bckg calculation and not checked again.
           int background_count=0;
           for (int j=0;j<nbins;j++)  // for one GTI selected (i), looping over all the others (j).
           {
                if (j!=i &&accepted_bin_bckg_vector[j]==1)
                {
                    background_count+= m_counts(pix_number, j);
                    alpha++;
                 }
           }

           background_bin_array[i] = background_count; //The background is averaged on the number of bins -1
           non = m_counts(pix_number, i); 
           noff = background_bin_array[i];
           alpha = (1./alpha);
          
           ///////////////////////////////////////////////////////////////////////////////// 
           // Compute sensitivity in Gaussian sigma
           double alpha1 = alpha + 1.0;
           double ntotal = non+noff;
           double arg1   = non/ntotal;
           double arg2   = noff/ntotal;
           double term1  = non * std::log((alpha1/alpha)*arg1);
           double term2  = noff * std::log(alpha1*arg2);
           sig  = sqrt(2.0 * (term1 + term2));
           /////////////////////////////////////////////////
           //sig_bin_vector[i]=sig;
           sig_histogram(i)=sig;
           excess_bin_vector[i]=non - alpha*noff;
           std::cout << "significance of the bin : " << i << " : " << sig << " - alpha : " << alpha << " - non: " << non << " - noff: " << noff << "- excess: " << non - alpha*noff <<  std::endl;
           if (sig>4.5 ) // if the bin is significant, it is removed from the bckg and we loop again.
           {
               accepted_bin_bckg_vector[i]=0;
               background_validated=false;
           }
        }

    }
}


/***********************************************************************//**
 * @brief Get the map index associated with a given time
 *
 * @param[in] time      Time
 * @return Index of map in counts cube
 ***************************************************************************/
int ctfindvar::time2inx(const GTime& time)
{
    // Set default return value
    int map_index = -1;

    // Loop over all GTIs
    for (int i=0; i<m_gti.size(); i++) {
        // Check if interval contains the time
        GGti interval = inx2gti(i);
        if (interval.contains(time)) {
            map_index = i;
            break;
        }
    }

    return map_index;
}


/***********************************************************************//**
 * @brief Return reference to significance cube object
 ***************************************************************************/
GGti ctfindvar::inx2gti(const int& index)
{
    if ((index < 0) || (index >= m_gti.size())) {
        // Index is invalid, so throw an error
        throw GException::invalid_value("ctfindvar::inx2gti(const int&)",
                                        "'index' parameter out of range");
    }

    return GGti(m_gti.tstart(index), m_gti.tstop(index));
}


/***********************************************************************//**
 * @brief Save something
 *
 * Saves something.
 ***************************************************************************/
void ctfindvar::save(void)
{
    log_header1(TERSE, "Saving results");

    // Filenames
    GFilename outcube((*this)["prefix"].string() + "countscube.fits");

    // Write counts cube
    log_value(TERSE, "Saving counts cube", outcube);
    m_counts.save(outcube, (*this)["clobber"].boolean());
    
    // Write individual source distributions

    // TODO: Your code goes here

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
void ctfindvar::init_members(void)
{
    // Initialise members
    m_counts.clear();
    m_gti.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] app search time variability tool.
 ***************************************************************************/
void ctfindvar::copy_members(const ctfindvar& app)
{
    // Copy attributes
    // TODO: Your code goes here

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void ctfindvar::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get application parameters
 *
 * @todo Implement method
 ***************************************************************************/
void ctfindvar::get_parameters(void)
{
    // Load the observations
    setup_observations(m_obs, false, true, false);

    // Create GTIs and counts cube
    init_gtis();
    init_cube();

    // Get the rest of the parameters
    if (read_ahead()) {
        (*this)["prefix"].string();
        (*this)["xsrc"].real();
        (*this)["ysrc"].real();
    }

    // Write parameters into logger
    log_parameters(TERSE);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill the cube from the events in the map
 *
 * Sets up the required event cube and time interval storage objects 
 * associated with this object
 ***************************************************************************/
void ctfindvar::fill_cube(GCTAObservation* obs)
{
    // Make sure that the observation holds a CTA event list. If this
    // is not the case then throw an exception.
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>
                                  (obs->events());
    if (events == NULL) {
        std::string msg = "CTA Observation does not contain an event "
                          "list. An event list is needed to fill the "
                          "counts cube.";
        throw GException::invalid_value(G_FILL_CUBE, msg);
    }

    // Get the RoI
    const GCTARoi& roi = events->roi();

    // Check for RoI sanity
    if (!roi.is_valid()) {
        std::string msg = "No valid RoI found in input observation "
                          "\""+obs->name()+"\". Run ctselect to specify "
                          "an RoI for this observation before running "
                          "ctbin.";
        throw GException::invalid_value(G_FILL_CUBE, msg);
    }

    // Initialise binning statistics
    int num_outside_roi  = 0;
    int num_invalid_wcs  = 0;
    int num_outside_map  = 0;
    int num_outside_time = 0;
    int num_in_map       = 0;

    // Fill counts sky map
    for (int i = 0; i < events->size(); ++i) {

        // Get event
        const GCTAEventAtom* event = (*events)[i];

        // Determine event time
        GTime evnt_time = event->time();
        int time_indx = time2inx(evnt_time);
        
        // Check if event is in the GTIs
        if ((time_indx > -1) && (time_indx < m_counts.nmaps())) {

            // Get event direction
            GSkyDir evnt_dir = GCTAInstDir(event->dir()).dir();

            // Skip event if it is outside the RoI
            if (roi.centre().dir().dist_deg(evnt_dir) > roi.radius()) {
                num_outside_roi++;
                continue;
            }

            // Determine sky pixel
            GSkyPixel pixel;
            try {
                pixel = m_counts.dir2pix(evnt_dir);
            }
            catch (std::exception &e) {
                num_invalid_wcs++;
                continue;
            }

            // Skip event if corresponding counts cube pixel is outside the
            // counts cube map range
            if (pixel.x() < -0.5 || pixel.x() > (m_counts.nx()-0.5) ||
                pixel.y() < -0.5 || pixel.y() > (m_counts.ny()-0.5)) {
                num_outside_map++;
                continue;
            }

            // Fill event in skymap
            #pragma omp critical(ctfindvar_fill_cube)
            m_counts(pixel, time_indx) += 1.0;

            // Increment number of maps
            num_in_map++;

        } else {
            // Event falls outside the events
            num_outside_time++;
            continue;
        }

    } // endfor: looped over all events

    // Log filling results
    #pragma omp critical(ctfindvar_fill_cube)
    {
        log_header3(TERSE, get_obs_header(obs));
        log_value(NORMAL, "Events in list", obs->events()->size());
        log_value(NORMAL, "Events in cube", num_in_map);
        log_value(NORMAL, "Events outside RoI", num_outside_roi);
        log_value(NORMAL, "Events with invalid WCS", num_invalid_wcs);
        log_value(NORMAL, "Events outside cube area", num_outside_map);
        log_value(NORMAL, "Events outside time bins", num_outside_time);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialize counts cube
 *
 * Sets up the required event cube and time interval storage objects 
 * associated with this object
 ***************************************************************************/
void ctfindvar::init_cube(void)
{
    // Create the basic skymap
    m_counts = create_map(m_obs);

    // Resize to the appropriate number of time intervals
    m_counts.nmaps(m_gti.size());
}


/***********************************************************************//**
 * @brief Initialize counts cube
 *
 * Sets up the required event cube and time interval storage objects 
 * associated with this object
 ***************************************************************************/
void ctfindvar::init_gtis(void)
{
    double tinterval = (*this)["tinterval"].real();
    GTime tstart("JD 999999999");     // large unreasonable Julian date
    GTime tstop("JD 0");              // small unreasonable Julian date

    // Set the start time
    if ((*this)["tmin"].is_valid()) {
        log_string(NORMAL, "Setting start time from command line");
        GTime tstart = (*this)["tmin"].time();
    } else {
        log_string(NORMAL, "Setting start time from observations");
        for (int o=0; o<m_obs.size(); o++) {
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[o]);
            
            // Check that the start time of the run is less than tstart
            if (obs->gti().tstart() < tstart) {
                tstart = obs->gti().tstart();
            }
        }
    }

    // Set the stop time
    if ((*this)["tmax"].is_valid()) {
        log_string(NORMAL, "Setting stop time from command line");
        GTime tstop = (*this)["tmax"].time();
    } else {

        log_string(NORMAL, "Setting stop time from observations");
        for (int o=0; o<m_obs.size(); o++) {
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[o]);
            
            // Check that the stop time of the run is larger than tstop
            if (obs->gti().tstop() > tstop) {
                tstop = obs->gti().tstop();
            }
        }
    }

    // Fill the number of good time intervals
    double tstart_sec = tstart.secs();
    double tstop_sec  = tstop.secs();
    int    bins       = (tstop_sec - tstart_sec) / tinterval + 0.5;

    // Log information
    log_value(NORMAL, "treference (mjd)", tstart.reference().mjdref());
    log_value(NORMAL, "tstart (sec)",  tstart.secs());
    log_value(NORMAL, "tstop  (sec)",  tstop.secs());
    log_value(NORMAL, "Total time",    tstop_sec-tstart_sec);
    log_value(NORMAL, "Time interval", tinterval);
    log_value(NORMAL, "nbins",         bins);

    // Throw if tstart == tstop
    if (tstart == tstop) {
        throw GException::invalid_value(G_INIT_GTIS,
                                        "Start time is equal to stop time");
    } else if (bins < 2) {
        std::stringstream msg;
        msg << "Method requires at least two time bins (" 
            << bins << " bins found).";
        throw GException::invalid_value(G_INIT_GTIS, msg.str());
    }

    for (int i=0; i<bins; i++) {
        
        // Get the next time interval
        GTime next_tstart(tstart + i*tinterval);
        GTime next_tstop(tstart + (i+1)*tinterval);
        GGti  next_gti(next_tstart, next_tstop);

        // Make sure there's an observation that overlaps with this GTI
        bool has_obs = false;
        for (int o=0; o<m_obs.size(); o++) {
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[o]);
            
            if (gtis_overlap(obs->gti(), next_gti)) {
                m_gti.extend(next_gti);
                break;
            }
        }
    }
}


/***********************************************************************//**
 * @brief Check if two GTIs overlap
 *
 * @param[in] gti1      First GTI
 * @param[in] gti2      Second GTI
 * @return Whether gti1 and gti2 overlap in time
 ***************************************************************************/
bool ctfindvar::gtis_overlap(const GGti& gti1, const GGti& gti2)
{
    bool overlap = false;
    // Loop over all GTIs
    if ((gti1.tstart() < gti2.tstart()) && (gti1.tstop() > gti2.tstart())) {
        overlap = true;
    } else if ((gti1.tstart() < gti2.tstop()) && (gti1.tstop() > gti2.tstop())) {
        overlap = true;
    } else if ((gti1.tstart() > gti2.tstart()) && (gti1.tstart() < gti2.tstop())) {
        overlap = true;
    }

    return overlap;
}
