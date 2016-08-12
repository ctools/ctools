/***************************************************************************
 *                        ctbin - Event binning tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @file ctbin.hpp
 * @brief Event binning tool definition
 * @author Juergen Knoedlseder
 */

#ifndef CTBIN_HPP
#define CTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include "ctobservation.hpp"

/* __Definitions _________________________________________________________ */
#define CTBIN_NAME    "ctbin"
#define CTBIN_VERSION "1.2.0"


/***********************************************************************//**
 * @class ctbin
 *
 * @brief Event binning tool
 *
 * This class bins event list(s) into a single counts cube. The class can
 * operate on predefined observation containers, on individual event list
 * FITS files, and on observation definition XML files.
 *
 * If multiple event lists are specified in the observation container or the
 * XML definition file, the class will merge these events into a single
 * counts cube.
 *
 * Results are stored in an observation container that can be written to disk
 * in form of a single FITS file. On output, the observation container will
 * have merged the input event lists into a single observation.
 *
 * WARNING: Note that the pointing direction of the counts cube will be set
 * to the skymap centre used for the counts cube definition. If usepnt=yes
 * is used, the pointing direction will be extracted from the first
 * observation encountered in the list. Ultimately, pointing direction
 * information should not be removed from the counts cube and exposure and
 * PSF cubes should be used for response computation. This is however not
 * yet fully implemented.
 ***************************************************************************/
class ctbin : public ctobservation {
public:
    // Constructors and destructors
    ctbin(void);
    explicit ctbin(const GObservations& obs);
    ctbin(int argc, char *argv[]);
    ctbin(const ctbin& app);
    virtual ~ctbin(void);

    // Operators
    ctbin& operator=(const ctbin& app);

    // Methods
    void                 clear(void);
    void                 run(void);
    void                 save(void);
    void                 publish(const std::string& name = "");
    const GCTAEventCube& cube(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctbin& app);
    void free_members(void);
    void get_parameters(void);
    void fill_cube(GCTAObservation* obs);
    void set_weights(GCTAObservation* obs);
    void obs_cube(void);

    // User parameters
    GFilename     m_outcube;  //!< Output counts map file name
    bool          m_usepnt;   //!< Use pointing instead of xref/yref parameters
    bool          m_publish;  //!< Publish counts cube?
    GChatter      m_chatter;  //!< Chattiness

    // Protected members
    GSkyMap       m_counts;   //!< Event cube counts
    GSkyMap       m_weights;  //!< Event cube weights
    GEbounds      m_ebounds;  //!< Energy boundaries
    GGti          m_gti;      //!< Good time intervals
    GCTAEventCube m_cube;     //!< Events cube (for cube() method)
    double        m_ontime;   //!< Total ontime
    double        m_livetime; //!< Total livetime
};


/***********************************************************************//**
 * @brief Return event cube
 *
 * @return Reference to event cube
 *
 * Returns a reference to the event cube.
 ***************************************************************************/
inline
const GCTAEventCube& ctbin::cube(void) const
{
    return m_cube;
}

#endif /* CTBIN_HPP */
