/***************************************************************************
 *          ctphase - Append phase information to CTA events file          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2022 by Joshua Cardenzana                           *
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
 * @file ctphase.hpp
 * @brief Event phase computation tool interface definition
 * @author Joshua Cardenzana
 */

#ifndef CTPHASE_HPP
#define CTPHASE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTPHASE_NAME "ctphase"


/***********************************************************************//**
 * @class ctphase
 *
 * @brief Event phase computation tool
 *
 * Computes the phase of each event using a temporal phase curve model
 ***************************************************************************/
class ctphase : public ctobservation {

public:
    // Constructors and destructors
    ctphase(void);
    explicit ctphase(const GObservations& obs);
    ctphase(int argc, char *argv[]);
    ctphase(const ctphase& app);
    virtual ~ctphase(void);

    // Operators
    ctphase& operator=(const ctphase& app);
    
    // Methods
    void clear(void);
    void process(void);
    void save(void);
    void publish(const std::string& name = "");

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctphase& app);
    void free_members(void);
    void get_parameters(void);
    void phase_events(GCTAObservation* obs);

    // Protected members
    GModelTemporalPhaseCurve m_phase; //!< Phase curve model
};


#endif /* CTPHASE_HPP */
