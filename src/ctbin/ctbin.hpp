/***************************************************************************
 *                        ctbin - Event binning tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2022 by Juergen Knoedlseder                         *
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
#define CTBIN_NAME "ctbin"


/***********************************************************************//**
 * @class ctbin
 *
 * @brief Event binning tool
 *
 * This class bins event list(s) into counts cubes. The class can
 * operate on predefined observation containers, on individual event list
 * FITS files, and on observation definition XML files.
 *
 * If multiple event lists are specified in the observation container or the
 * XML definition file, the class will upon request stack these events into a
 * single counts cube.
 *
 * If the hidden parameter stack=no is used, one counts cube is generated for
 * each individual observation. These cubes are saved to a path that can be
 * specified by the hidden parameter prefix, followed by the corresponding
 * observation id. An observation definition XML file containing the paths
 * to the newly generated counts cubes is written to path given in the 
 * parameter outcube.
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
    void                 process(void);
    int                  cubes(void) const;
    const GCTAEventCube& cube(const int& index = 0) const;
    void                 save(void);
    void                 publish(const std::string& name = "");

protected:
    // Protected methods
    void    init_members(void);
    void    copy_members(const ctbin& app);
    void    free_members(void);
    void    get_parameters(void);
    void    init_sky_dir_cache(void);
    GSkyMap create_cube(const GCTAObservation* obs);
    void    fill_cube(const GCTAObservation* obs, GSkyMap& counts);
    void    set_weights(const GCTAObservation* obs, GSkyMap& weights);
    void    obs_cube_stacked(void);

    // User parameters
    bool        m_stack;    //!< Output one stacked cube or multiple cubes
    std::string m_prefix;   //!< Prefix for multiple cubes
    bool        m_publish;  //!< Publish counts cube?
    GChatter    m_chatter;  //!< Chattiness

    // Protected members
    std::vector<GCTAEventCube> m_cubes;    //!< Event cubes
    std::vector<GSkyMap>       m_counts;   //!< List of event cube counts
    std::vector<GSkyMap>       m_weights;  //!< List of event cube weights
    GSkyDir                    m_mean_pnt; //!< Mean pointing
    GEbounds                   m_ebounds;  //!< Energy boundaries
    GGti                       m_gti;      //!< Stacked Good time intervals
    double                     m_ontime;   //!< Total ontime
    double                     m_livetime; //!< Total livetime

    // Cache members
    std::vector<GSkyDir>       m_dirs;     //!< Cached GSkyDir for all pixels
};


/***********************************************************************//**
 * @brief Return number of event cubes
 *
 * @return Number of event cubes.
 *
 * Returns number of event cubes.
 ***************************************************************************/
inline
int ctbin::cubes(void) const
{
    return ((int)m_cubes.size());
}

#endif /* CTBIN_HPP */
