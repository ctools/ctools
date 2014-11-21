/***************************************************************************
 *                ctbutterfly - butterfly calculation tool                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Michael Mayer                                    *
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
 * @file ctbutterfly.hpp
 * @brief Butterfly calculation tool interface definition
 * @author Michael Mayer
 */

#ifndef CTBUTTERFLY_HPP
#define CTBUTTERFLY_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "ctool.hpp"

/* __Definitions _________________________________________________________ */
#define CTBUTTERFLY_NAME    "ctbutterfly"
#define CTBUTTERFLY_VERSION "00-02-00"


/***********************************************************************//**
 * @class ctbutterfly
 *
 * @brief Butterfly calculation tool
 *
 * This class computes the confidence interval of a fitted spectrum
 * (butterfly) for a given set of observations.
 *
 * The class operates on predefined observation containers, an individual
 * event list or an observation definition XML file.
 *
 * During the computation the covariance matrix of the fit is used to
 * propagate uncertainties and their correlations through the entire energy
 * range. The output is saved as an ascii files containing the confidence
 * band boundaries
 ***************************************************************************/
class ctbutterfly : public ctool {

public:
    // Constructors and destructors
    ctbutterfly(void);
    explicit ctbutterfly(const GObservations& obs);
    ctbutterfly(int argc, char *argv[]);
    ctbutterfly(const ctbutterfly& app);
    virtual ~ctbutterfly(void);

    // Operators
    ctbutterfly& operator=(const ctbutterfly& app);

    // Methods
    void                 clear(void);
    void                 run(void);
    void                 save(void);
    const GObservations& obs(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctbutterfly& app);
    void free_members(void);
    void get_parameters(void);

    // User parameters
    std::string m_srcname;    //!< Name of source to compute butterfly
    std::string m_outfile;    //!< Output ascii file
    GEbounds    m_ebounds;    //!< Energy binning definition

      // Protected members
    GObservations       m_obs;        //!< Observation container
    GMatrixSparse       m_covariance; //!< Covariance matrix
    std::vector<double> m_energies;   //!< Energy values for storage
    std::vector<double> m_fluxes;     //!< Flux values per energy bin
    std::vector<double> m_errors;     //!< Flux errors per energy bin
};

/***********************************************************************//**
 * @brief Return observation container
 *
 * @return Reference to observation container
 ***************************************************************************/
inline
const GObservations& ctbutterfly::obs(void) const
{
    return m_obs;
}

#endif /* CTBUTTERFLY_HPP */
