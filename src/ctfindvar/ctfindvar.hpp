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
    void get_parameters(void);
    void get_variability_sig(void)
    void get_variability_sig(const int& pix_number, const int& nbins, GNdarray& sig_histogram)

    // Protected members
    GSkyMap m_counts;                    //!< Counts for each time interval
    GGti    m_gti;                            //!< List of time intervals
};


/***********************************************************************//**
 * @brief Return reference to significance cube object
 ***************************************************************************/
inline
const GSkyMap& ctfindvar::counts(void)
{
    return m_counts;
}


#endif /* CTFINDVAR_HPP */
