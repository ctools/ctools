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

/* __ OpenMP section _____________________________________________________ */
//#ifdef _OPENMP
//#include <omp.h>
//#endif

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
 *
 * Constructs empty search time variability tool.
 ***************************************************************************/
ctfindvar::ctfindvar(void) : ctobservation(CTFINDVAR_NAME, VERSION)
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
ctfindvar::ctfindvar(const GObservations& obs) : ctobservation(CTFINDVAR_NAME, VERSION, obs)
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
ctfindvar::ctfindvar(int argc, char *argv[]) : ctobservation(CTFINDVAR_NAME, VERSION, argc, argv)
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
ctfindvar::ctfindvar(const ctfindvar& app) : ctobservation(app)
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

    // TODO: Your code goes here

   // //////////////////////////////////////////////////////////////////
   // double max_sig=0;
   // GNdarray GNdarray_temp(NBINS);
   // for (int i=0, i<number of pixels; i++)
   // {    
   //     get_variability_sig(pix_number, GNdarray_temp) 
   //     if (position == source position)
   //     {
   //         GNdarray_source= GNdarray_temp;
   //     }
   //     else
   //     {
   //         if (max(GNdarray_temp) > max_sig)
   //         {
   //             GNdarray_max = GNdarray_temp;
   //         }
   //     }
   // }
    //////////////////////////////////////////////////////////////////
    // Return
    return;
}
void ctfindvar::get_variability_sig(const int& pix_number, const int& nbins, GNdarray& sig_histogram)
{
    //const int nbins = hist->GetNbinsX();
    //const int nbins = VECTOR.size() 
    //const int nbins = GSKYMAP.maps(); 
    int background_bin_array[nbins];
    bool background_validated=false;
    int non, noff;
    double alpha, sig;

    for (int i=0;i<nbins;i++)
    {
        background_bin_array[i]=0;
        accepted_bin_bckg_vector.push_back(1);
        //sig_bin_vector.push_back(0);
        excess_bin_vector.push_back(0);
    }

    while (background_validated==false)
    {
        background_validated=true;

        for (int i=0; i< nbins; i++) //looping over all the bins of the histo
        {
           alpha=0;
           if (accepted_bin_bckg_vector[i]==0) continue;     //the run is discared from bck calculation and not checked again.
           int background_count=0;
           for (int j=0;j<nbins;j++)  // for one bin selected (i), looping over all the others (j).
           {
                if (j!=i &&accepted_bin_bckg_vector[j]==1)
                {
                    //background_count+=hist->GetBinContent(j+1);
                    //background_count+=VECTOR(j+1).first(pix_number);
                    background_count+=GSKYMAP(pix_number, j);
                    alpha++;
                 }
           }

           background_bin_array[i] = background_count; //The background is averaged on the number of bins -1
           //non = hist->GetBinContent(i+1);
           //non = VECTOR(i+1).first(pix_number); 
           non = GSKYMAP(pix_number, i); 
           noff = background_bin_array[i];
           alpha = (1./alpha);
           sig = Utilities::Statistics::Significance(non,noff,alpha);
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
 * @brief Save something
 *
 * Saves something.
 ***************************************************************************/
void ctfindvar::save(void)
{
    // Write header
    log_header1(TERSE, "Save something");

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


    // Write parameters into logger
    log_parameters(TERSE);

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
    GTime tstart("JD 999999999");          // large unreasonable Julian date
    GTime tstop("JD 0");                   // small unreasonable Julian date

    // Set the start time
    if ((*this)["tmin"].is_valid()) {
        GTime tstart = (*this)["tmin"].time();
    } else {
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
        GTime tstop = (*this)["tmax"].time();
    } else {
        for (int o=0; o<m_obs.size(); o++) {
            GCTAObservation* obs = dynamic_cast<GCTAObservation*>(m_obs[o]);
            
            // Check that the start time of the run is less than tstart
            if (obs->gti().tstop() < tstop) {
                tstop = obs->gti().tstop();
            }
        }
    }

    // ToDO: This is wrong only making 1 bin
    m_gti = GGti(tstart, tstop);
}