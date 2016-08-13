/***************************************************************************
 *             ctobservation - Base class for observation tools            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file ctobservation.hpp
 * @brief Observation tool base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef CTOBSERVATION_HPP
#define CTOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "ctool.hpp"

/* __Definitions _________________________________________________________ */


/***********************************************************************//**
 * @class ctobservation
 *
 * @brief Base class for observation tools
 *
 * This is the baseclass for observation tools. Observation tools are ctools
 * that hold an observation container.
 ***************************************************************************/
class ctobservation : public ctool {

public:
    // Constructors and destructors
    ctobservation(const std::string& name, const std::string& version);
    ctobservation(const std::string& name, const std::string& version,
                  const GObservations& obs);
    ctobservation(const std::string& name, const std::string& version,
                  int argc, char* argv[]);
    ctobservation(const ctobservation& app);
    virtual ~ctobservation(void);

    // Operators
    ctobservation& operator=(const ctobservation& app);

    // Pure virtual methods
    virtual void clear(void) = 0;
    virtual void run(void) = 0;
    virtual void save(void) = 0;

    // Methods
    const GObservations& obs(void) const;

#ifndef SWIG
protected:
#endif
    // Iterator methods and members
    GCTAObservation*       first_unbinned_observation(void);
    GCTAObservation*       next_unbinned_observation(void);
    const GCTAObservation* first_unbinned_observation(void) const;
    const GCTAObservation* next_unbinned_observation(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctobservation& app);
    void free_members(void);

    // Protected members
    GObservations m_obs;            //!< Observation container

private:
    // Private members
    mutable int m_index_unbinned;   //!< Current index of unbinned observation
};


/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container.
 *
 * Returns a reference to the observation container.
 ***************************************************************************/
inline
const GObservations& ctobservation::obs(void) const
{
    return m_obs;
}


/***********************************************************************//**
 * @brief Return first unbinned CTA observation
 *
 * @return Pointer to first unbinned CTA observation
 *
 * Returns a pointer to the first unbinned CTA observation in the container.
 * If no CTA observation exists a NULL pointer is returned.
 *
 * The method calls next_unbinned_observation(). See the method for details.
 ***************************************************************************/
inline
GCTAObservation* ctobservation::first_unbinned_observation(void)
{
    return (const_cast<GCTAObservation*>
            (static_cast<const ctobservation&>
             (*this).first_unbinned_observation()));
}


/***********************************************************************//**
 * @brief Return next unbinned CTA observation
 *
 * @return Pointer to next unbinned CTA observation
 *
 * Returns a pointer to the next unbinned CTA observation in the container.
 * If no CTA observation exists any more a NULL pointer is returned.
 *
 * The method writes for each encountered observation a level 3 header into
 * the logger. It will also signal when an observation was skipped because
 * it either was not a CTA observation or not an unbinned observation.
 ***************************************************************************/
inline
GCTAObservation* ctobservation::next_unbinned_observation(void)
{
    return (const_cast<GCTAObservation*>
            (static_cast<const ctobservation&>
             (*this).next_unbinned_observation()));
}

#endif /* CTOBSERVATION_HPP */
