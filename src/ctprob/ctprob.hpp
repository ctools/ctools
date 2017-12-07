/***************************************************************************
 *          ctprob - Computes event probability for a given model          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Leonardo Di Venere                               *
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
 * @file ctprob.hpp
 * @brief Event probability computation tool interface definition
 * @author Leonardo Di Venere
 */

#ifndef CTPROB_HPP
#define CTPROB_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTPROB_NAME "ctprob"


/***********************************************************************//**
 * @class ctprob
 *
 * @brief Event probability computation tool
 ***************************************************************************/
class ctprob : public ctobservation {

public:
    // Constructors and destructors
    ctprob(void);
    explicit ctprob(const GObservations& obs);
    ctprob(int argc, char *argv[]);
    ctprob(const ctprob& app);
    virtual ~ctprob(void);

    // Operators
    ctprob& operator=(const ctprob& app);

    // Methods
    void clear(void);
    void run(void);
    void save(void);
    void publish(const std::string& name = "");

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctprob& app);
    void free_members(void);
    void get_parameters(void);
    void evaluate_probability(GCTAObservation* obs);

    // User parameters
    bool     m_apply_edisp; //!< Apply energy dispersion?
    bool     m_publish;     //!< Publish event list?
    GChatter m_chatter;     //!< Chattiness
};

#endif /* CTPROB_HPP */
