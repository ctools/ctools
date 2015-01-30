/***************************************************************************
 *                  ctobssim - Observation simulator tool                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2015 by Juergen Knoedlseder                         *
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
#include "ctool.hpp"

/* __Definitions _________________________________________________________ */
#define CTOBSSIM_NAME    "ctobssim"
#define CTOBSSIM_VERSION "00-10-00"


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
class ctobssim : public ctool {

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
    void                 clear(void);
    void                 run(void);
    void                 save(void);
    const GObservations& obs(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctobssim& app);
    void free_members(void);
    void get_parameters(void);
    void simulate_source(GCTAObservation* obs,
                         const GModels&   models,
                         GRan&            ran, 
                         GLog*            wrklog = NULL);
    void simulate_background(GCTAObservation* obs,
                             const GModels&   models,
                             GRan&            ran,
                             GLog*            wrklog = NULL);
    void save_fits(void);
    void save_xml(void);

    // User parameters
    std::string m_outevents;   //!< Output events file
    std::string m_prefix;      //!< Prefix for multiple event lists
    int         m_seed;        //!< Random number generator seed
    bool        m_apply_edisp; //!< Apply energy dispersion?

    // Protected members
    mutable GObservations m_obs;              //!< Observation container
    mutable bool          m_save_and_dispose; //!< Save and dispose immediately
    double                m_area;             //!< Surface area for simulation (cm2)
    int                   m_max_photons;      //!< Maximum number of photons/slice
    std::vector<GRan>     m_rans;             //!< Random number generators
    int                   m_event_id;         //!< Event identifier
};


/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GObservations& ctobssim::obs(void) const
{
    return m_obs;
}

#endif /* CTOBSSIM_HPP */
