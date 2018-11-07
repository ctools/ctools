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

    // creating GSkyDir to get the position of the source and each pixels

    // Extract which pixel in the map each source is located in
    std::vector<int> srcInxPix;
    if (m_inmodel.size() > 0) {
        for (int mod=0; mod<m_inmodel.size(); mod++) {

            // Extract the model
            GModelSky* model = dynamic_cast<GModelSky*>(m_inmodel[mod]);
            if (model != NULL) {
                // Get the source position
                GModelSpatial* model_spatial = model->spatial();
                GSkyRegion*    model_region  = model_spatial->region();
                GSkyDir        srcSkyDir     = 
                        dynamic_cast<GSkyRegionCircle*>(model_region)->centre();

                // Store the map index of this source
                int src_map_index = m_counts.dir2inx(srcSkyDir);
                srcInxPix.push_back(src_map_index);
            } else {
                m_inmodel.remove(mod);
                mod--;
            }
        }
    } 
    // ... otherwise only one source position is queried
    else {
    GSkyDir srcSkyDir;
    if ((*this)["coordsys"].string()=="CEL")
    {
        srcSkyDir.radec_deg((*this)["xsrc"].real(),(*this)["ysrc"].real());
    }
    else
    {
        srcSkyDir.lb_deg((*this)["xsrc"].real(),(*this)["ysrc"].real());
    }

        // Store the index of this source
        int src_map_index = m_counts.dir2inx(srcSkyDir);
        srcInxPix.push_back(src_map_index);
    }

    //preparing the histogram with significance evolution for the source center and the highest sig-pix
    GNdarray pixSig(nbins), pixSigMax(nbins); //Storing the significance for each GTI of the pixels
    m_pixsigsrc = GNdarray(srcInxPix.size(), nbins);
    m_pixsigmax = GNdarray(1,nbins);

    //Prepare the final skymap with the max significance of each pixel
    double  max_sig=0;
    GSkyDir max_sig_dir;

    //looping over all the pixels in the cube
    for (int pix_number=0; pix_number<m_counts.npix(); pix_number++)
    {    
        double total_counts=0;
        for (int k=0;k<m_counts.nmaps();k++)
        {
            total_counts += m_counts(pix_number, k);
        }
        #ifdef G_DEBUG
        if(pix_number%20==0)
        {
            std::cout << "Pixel number " << pix_number << " has a total number of counts of: " << total_counts << std::endl;
        }
        #endif

        //Getting the variability significance for the pixel
        get_variability_sig(pix_number,nbins, pixSig);

        //Getting the significance evolution for the source
        for (int src=0; src<srcInxPix.size(); src++) {

            // Store the distribution if the source is located at this position
            if (srcInxPix[src] == pix_number) {
            for (int i=0; i<nbins; i++) {

                    m_pixsigsrc(src,i) = pixSig(i);
            }
            #ifdef G_DEBUG
            std::cout << "checking pixel number of the source of interest" << pix_number << std::endl;
            std::cout << "number of counts in pixel of the source of interest" << total_counts << std::endl;
            pixSigSrc=pixSig;
            std::cout << "checking pixel number " << pix_number << " with total number of counts of: " << total_counts << std::endl;
            std::cin.ignore();
            std::cout << "checking pixel number of the source of interest: " << pix_number << " - with total number of counts of: " << total_counts << std::endl;
            #endif
        }
        }
        
            // Getting the evolution for the pix with highest significance
            if (max(pixSig) > max_sig)  
            {
            max_sig = max(pixSig);
            // Store information for pixel with the maximum significance
            for (int i=0; i<nbins; i++) {
                m_pixsigmax(0,i) = pixSig(i);
            }

            // Update the pixel position with the maximum sigma
            max_sig     = max(pixSig);
            max_sig_dir = m_counts.inx2dir(pix_number);
        }

        //storing sig in skymap
        pixVarSig(pix_number) = max(pixSig);
    }

    log_header1(NORMAL, "Analysis finished");
    log_value(NORMAL, "Maximum sigma", max_sig);
    log_value(NORMAL, "Max pixel RA", max_sig_dir.ra_deg());
    log_value(NORMAL, "Max pixel DEC", max_sig_dir.dec_deg());

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
    double non, noff;
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
           if (noff == 0.0) {
               sig = 0.0;
           } else if (non > 0.0) {
               double term1  = non * std::log((alpha1/alpha)*arg1);
               double term2  = noff * std::log(alpha1*arg2);
               sig  = std::sqrt(2.0 * (term1 + term2));
           } else {
               sig = std::sqrt(2.0 * noff * std::log(alpha1));
           }

           // Specify the sign of the significance
           sig *= (non < alpha*noff) ? -1.0 : 1.0;

           /////////////////////////////////////////////////
           //sig_bin_vector[i]=sig;
           sig_histogram(i)=sig;
           excess_bin_vector[i]=non - alpha*noff;
           #ifdef G_DEBUG
           //std::cout << "significance of the bin : " << i << " : " << sig << " - alpha : " << alpha << " - non: " << non << " - noff: " << noff << "- excess: " << non - alpha*noff <<  std::endl;
           #endif
           if (sig>4.5 ) // if the bin is significant, it is removed from the bckg and we loop again.
           {
               accepted_bin_bckg_vector[i]=0;
               background_validated=false;
           }
        }

    }
    if (pix_number==19900)
    {
        for (int k=0;k<sig_histogram.size();k++)
        {
            std::cout << "value of the sig of the " << k << "th slice =" << sig_histogram(k) << std::endl;
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
    std::string prefix((*this)["prefix"].string());
    GFilename outcube(prefix + "countscube.fits");
    GFilename peaksigmap(prefix + "peaksigmap.fits");

    // Write counts cube
    log_value(TERSE, "Saving counts cube", outcube);
    m_counts.save(outcube, (*this)["clobber"].boolean());
    
    // Write the most significant values for each pixel
    log_value(TERSE, "Saving maximum significances map", peaksigmap);
    m_peaksigmap.save(peaksigmap, (*this)["clobber"].boolean());

    // Write individual source distributions
    write_srchist();

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
    m_inmodel.clear();
    m_peaksigmap.clear();
    m_pixsigsrc.clear();
    m_tstart.clear();
    m_tstop.clear();

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

        if ((*this)["inmodel"].is_valid()) {
            m_inmodel = GModels((*this)["inmodel"].filename());
        } else {
            // Get the desired source position from user
            (*this)["xsrc"].real();
            (*this)["ysrc"].real();

            // Get the file format for the output source histograms
            (*this)["histtype"].string();
        }
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
            // Event falls outside the specified time range
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
    m_tstart = GTime("JD 999999999");     // large unreasonable Julian date
    m_tstop  = GTime("JD 0");              // small unreasonable Julian date

    // Set the start time
    if ((*this)["tmin"].is_valid()) {
        log_string(NORMAL, "Setting start time from command line");
        m_tstart = (*this)["tmin"].time();
    } else {
        log_string(NORMAL, "Setting start time from observations");
        for (int o=0; o<m_obs.size(); o++) {
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[o]);
            
            // Check that the start time of the run is less than tstart
            if (obs->gti().tstart() < m_tstart) {
                m_tstart = obs->gti().tstart();
            }
        }
    }

    // Set the stop time
    if ((*this)["tmax"].is_valid()) {
        log_string(NORMAL, "Setting stop time from command line");
        m_tstop = (*this)["tmax"].time();
    } else {

        log_string(NORMAL, "Setting stop time from observations");
        for (int o=0; o<m_obs.size(); o++) {
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[o]);
            
            // Check that the stop time of the run is larger than tstop
            if (obs->gti().tstop() > m_tstop) {
                m_tstop = obs->gti().tstop();
            }
        }
    }

    // Fill the number of good time intervals
    double tstart_sec = m_tstart.secs();
    double tstop_sec  = m_tstop.secs();
    int    bins       = (tstop_sec - tstart_sec) / tinterval + 0.5;

    // Log information
    log_value(NORMAL, "treference (mjd)", m_tstart.reference().mjdref());
    log_value(NORMAL, "tstart (sec)",  tstart_sec);
    log_value(NORMAL, "tstop  (sec)",  tstop_sec);
    log_value(NORMAL, "Total time",    tstop_sec-tstart_sec);
    log_value(NORMAL, "Time interval", tinterval);
    log_value(NORMAL, "nbins",         bins);

    // Throw if tstart == tstop
    if (m_tstart == m_tstop) {
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
        GTime next_tstart(m_tstart + i*tinterval);
        GTime next_tstop(m_tstart + (i+1)*tinterval);
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


/***********************************************************************//**
 * @brief Write individual histograms in CSV format
 ***************************************************************************/
void ctfindvar::write_srchist(void)
{
    // Setup data storage for time and source significances
    double   tinterval = (*this)["tinterval"].real();
    int      nentries = (m_tstop-m_tstart) / tinterval + 0.5;
    int      nsources = m_pixsigsrc.shape()[0];
    GNdarray time_info(2,nentries);
    GNdarray src_info(nsources, nentries);

    // Store the full set of information
    for (int i=0; i<nentries; i++) {
        GTime tstart_bin(m_tstart + i * tinterval);
        GTime tstop_bin(m_tstart + (i+1) * tinterval);
        GTime midpoint((tstop_bin.secs()+tstart_bin.secs())*0.5);

        // Store the time of this bin
        time_info(0,i) = tstart_bin.mjd();
        time_info(1,i) = tstop_bin.mjd();

        // Get the index associated with this time
        int index = time2inx(midpoint);

        // Loop over every source
        for (int src=0; src<nsources; src++) {
        if (index >= 0) {
                src_info(src,i) = m_pixsigsrc(src,index);
        } else {
                src_info(src,i) = 0.0;
            }
        }
    }

    // Call the appropriate method for writing the data
    if ((*this)["histtype"].string() == "CSV") {
        write_srchist_csv(time_info, src_info);
    } else {
        write_srchist_fits(time_info, src_info);
    }
}

/***********************************************************************//**
 * @brief Write individual histograms in CSV format
 ***************************************************************************/
void ctfindvar::write_srchist_csv(const GNdarray& time_info,
                                  const GNdarray& src_info)
{
    log_string(NORMAL,"Storing histograms in CSV format is not currently supported.");
}


/***********************************************************************//**
 * @brief Write individual histograms in FITS format
 * 
 * @param[in] time_info     Start and stop time in MJD for each bin
 * @param[in] src_info      Significance distribution for each source

 * This method stores the extracted significance for each source at each time
 * interval between tstart and tstop. Each time bin contains information on
 * the start and stop time of that bin (in MJD) and the extracted variability
 * significance of the source in that time bin.
 ***************************************************************************/
void ctfindvar::write_srchist_fits(const GNdarray& time_info,
                                   const GNdarray& src_info)
{
    log_string(NORMAL, "Storing histograms in FITS format");

    // Create the table
    int           nrows = time_info.shape()[1];
    GFitsBinTable table(nrows);
    table.extname("SRCSIGNIF");

    // Append the start and stop time columns
    GFitsTableFloatCol start_mjd("TSTART", nrows);
    GFitsTableFloatCol stop_mjd("TSTOP", nrows);

    for (int i=0; i<nrows; i++) {
        start_mjd(i) = time_info(0,i);
        stop_mjd(i)  = time_info(1,i);
    }
    table.append(start_mjd);
    table.append(stop_mjd);

    // Append the distribution for each source
    int src_count = src_info.shape()[0];
    for (int src=0; src<src_count; src++) {
        std::string src_name("SOURCE");
        if (m_inmodel.size() > 0) {
            src_name = m_inmodel[src]->name();
        }
        
        GFitsTableFloatCol src_sig(src_name, nrows);

        for (int i=0; i<nrows; i++) {
            src_sig(i) = src_info(src, i);
        }
        table.append(src_sig);
    }

    // Write the table to a file
    GFilename outfile((*this)["prefix"].string() + "srcsig.fits");
    GFits     fitsfile;
    fitsfile.append(table);
    fitsfile.saveto(outfile, (*this)["clobber"].boolean() );
}
