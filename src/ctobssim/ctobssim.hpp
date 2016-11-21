/***************************************************************************
 *                  ctobssim - Observation simulator tool                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2016 by Juergen Knoedlseder                         *
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
 * @file ctobssim.hpp
 * @brief Observation simulator tool definition
 * @author Juergen Knoedlseder
 */

#ifndef CTOBSSIM_HPP
#define CTOBSSIM_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTOBSSIM_NAME    "ctobssim"
#define CTOBSSIM_VERSION "1.2.0"


/***********************************************************************//**
 * @class ctobssim
 *
 * @brief Observation simulator tool
 *
 * This class simulates CTA observation(s) using Monte Carlo sampling of the
 * source and background models. The class supports simulation of data of
 * multiple CTA observations in one shot. If multiple CTA observations are
 * processed and the save method is called, events FITS files will be written
 * for each observation.
 ***************************************************************************/
class ctobssim : public ctobservation {

public:
    // Constructors and destructors
    ctobssim(void);
    explicit ctobssim(const GObservations& obs);
    ctobssim(int argc, char *argv[]);
    ctobssim(const ctobssim& app);
    virtual ~ctobssim(void);

    // Operators
    ctobssim& operator=(const ctobssim& app);

    // Methods
    void          clear(void);
    void          run(void);
    void          save(void);
    void          publish(const std::string& name = "");
    const double& max_rate(void) const;
    void          max_rate(const double& max_rate);

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const ctobssim& app);
    void        free_members(void);
    void        get_parameters(void);
    void        simulate_source(GCTAObservation* obs,
                                const GModels&   models,
                                GRan&            ran,
                                GLog*            wrklog = NULL);
    void        simulate_interval(GCTAObservation*       obs,
                                  const GCTAResponseIrf* rsp,
                                  GCTAEventList*         events,
                                  const GModels&         models,
                                  const GTime&           tmin,
                                  const GTime&           tmax,
                                  const GEnergy&         etrue_min,
                                  const GEnergy&         etrue_max,
                                  const GEnergy&         ereco_min,
                                  const GEnergy&         ereco_max,
                                  const GSkyDir&         dir,
                                  const double&          rad,
                                  const double&          area,
                                  GRan&                  ran,
                                  GLog*                  wrklog,
                                  int&                   indent,
                                  std::vector<int>&      nphotons,
                                  std::vector<int>&      nevents);
    void        simulate_time_slice(GCTAObservation*       obs,
                                    const GCTAResponseIrf* rsp,
                                    GCTAEventList*         events,
                                    const GModelSky*       model,
                                    const GTime&           tstart,
                                    const GTime&           tstop,
                                    const GEnergy&         etrue_min,
                                    const GEnergy&         etrue_max,
                                    const GEnergy&         ereco_min,
                                    const GEnergy&         ereco_max,
                                    const GSkyDir&         dir,
                                    const double&          rad,
                                    const double&          area,
                                    GRan&                  ran,
                                    GLog*                  wrklog,
                                    int&                   indent,
                                    int&                   nphotons,
                                    int&                   nevents);
    GEbounds    get_ebounds(const GEbounds& ebounds) const;
    double      get_area(GCTAObservation* obs,
                         const GEnergy&   emin,
                         const GEnergy&   emax) const;
    double      get_model_flux(const GModelSky* model,
                               const GEnergy&   emin,
                               const GEnergy&   emax,
                               const GSkyDir&   centre,
                               const double&    radius,
                               const int&       indent,
                               GLog*            wrklog);
    void        simulate_background(GCTAObservation* obs,
                                    const GModels&   models,
                                    GRan&            ran,
                                    GLog*            wrklog = NULL);
    void        save_fits(void);
    void        save_xml(void);
    std::string outfile(const int& index);

    // User parameters
    std::string m_outevents;   //!< Output events file
    std::string m_prefix;      //!< Prefix for multiple event lists
    int         m_startindex;  //!< Start index for multiple event lists
    int         m_seed;        //!< Random number generator seed
    int         m_eslices;     //!< Number of energy slices
    bool        m_apply_edisp; //!< Apply energy dispersion?

    // Protected members
    mutable bool          m_save_and_dispose; //!< Save and dispose immediately
    int                   m_max_photons;      //!< Maximum number of photons/slice
    double                m_max_rate;         //!< Maximum photon rate
    std::vector<GRan>     m_rans;             //!< Random number generators
    int                   m_event_id;         //!< Event identifier
};


/***********************************************************************//**
 * @brief Return maximum photon rate
 *
 * @return Reference to maximum photon rate.
 ***************************************************************************/
inline
const double& ctobssim::max_rate(void) const
{
    return (m_max_rate);
}


/***********************************************************************//**
 * @brief Set maximum photon rate
 *
 * @param[in] max_rate Maximum photon rate.
 ***************************************************************************/
inline
void ctobssim::max_rate(const double& max_rate)
{
    m_max_rate = max_rate;
    return;
}

#endif /* CTOBSSIM_HPP */
