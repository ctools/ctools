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
 * @file ctfindvar.hpp
 * @brief search time variability tool definition
 * @author Simon Bonnefoy
 */

#ifndef CTFINDVAR_HPP
#define CTFINDVAR_HPP

/* __ Includes ___________________________________________________________ */
#include "ctobservation.hpp"
#include "GSkyMap.hpp"
#include "GGti.hpp"


/* __Definitions _________________________________________________________ */
#define CTFINDVAR_NAME "ctfindvar"


/***********************************************************************//**
 * @class ctfindvar
 *
 * @brief search time variability tool
 *
 * @todo Add tool description.
 ***************************************************************************/
class ctfindvar : public ctobservation {
public:
    // Constructors and destructors
    ctfindvar(void);
    explicit ctfindvar(const GObservations& obs);
    ctfindvar(int argc, char *argv[]);
    ctfindvar(const ctfindvar& app);
    virtual ~ctfindvar(void);

    // Operators
    ctfindvar& operator=(const ctfindvar& app);

    // Methods
    void clear(void);
    void run(void);
    void save(void);

    // Get the information on the time interval from 
    int            time2inx(const GTime& time);
    GGti           inx2gti(const int& indx);
    const GSkyMap& counts(void);

protected:
    // Protected methods
    void init_cube(void);
    void init_gtis(void);
    void init_members(void);
    void copy_members(const ctfindvar& app);
    void free_members(void);
    void fill_cube(GCTAObservation* obs);
    void fill_alpha_vector(const int&           pix_number,
                           std::vector<double>& alpha_vector);
    void get_parameters(void);
    void get_variability_sig(const int& pix_number, const int& nbins, GNdarray& sig_histogram);
    double gti_overlap(const GGti& gti1, const GGti& gti2);

    // Writing methods for individual source histograms
    void write_srchist(void);
    void write_srchist_csv(const GNdarray& time_info,
                           const GNdarray& max_pixel_info,
                           const GNdarray& src_info);
    void write_srchist_fits(const GNdarray& time_info,
                            const GNdarray& max_pixel_info,
                            const GNdarray& src_info);

    // Protected members
    GSkyMap           m_counts;        //!< Counts for each time interval
    std::vector<GGti> m_gti;           //!< List of time intervals
    GModels           m_inmodel;       //!< List of models for source positions
    GSkyDir           m_max_sig_dir;   //!< Sky direction associated with maximum significance
    double            m_minoff;        //!< Minimum counts for use in significance calculation
    double            m_sig_threshold; //!< Minimum significance required to set source as variable 
    GSkyMap           m_peaksigmap;    //!< Skymap holding the maximum significance
    GNdarray          m_pixsigsrc;     //!< Store distributions of the source significances
    GNdarray          m_pixsigmax;     //!< Store distribution for pixel with max significance
    GTime             m_tstart;        //!< Start time for variability study
    GTime             m_tstop;         //!< Stop time for variability study
    GEnergy           m_emin;          //!< Minimum energy for events
    GEnergy           m_emax;          //!< Maximum energy for events

};


/***********************************************************************//**
 * @brief Return reference to significance cube object
 * 
 * @return Reference to counts cube
 ***************************************************************************/
inline
const GSkyMap& ctfindvar::counts(void)
{
    return m_counts;
}


#endif /* CTFINDVAR_HPP */
