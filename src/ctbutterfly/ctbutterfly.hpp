/***************************************************************************
 *                ctbutterfly - butterfly calculation tool                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Michael Mayer                               *
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
#include "ctlikelihood.hpp"

/* __Definitions _________________________________________________________ */
#define CTBUTTERFLY_NAME    "ctbutterfly"
#define CTBUTTERFLY_VERSION "1.1.0"


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
class ctbutterfly : public ctlikelihood {

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
    void clear(void);
    void run(void);
    void save(void);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctbutterfly& app);
    void free_members(void);
    void get_parameters(void);
    void check_model(void);
    void eigenvectors(const double& a,
                      const double& b,
                      const double& c,
                      const double& d,
                      double*       lambda1,
                      double*       lambda2,
                      GVector*      vector1,
                      GVector*      vector2);

    // User parameters
    std::string m_srcname;      //!< Name of source to compute butterfly
    double      m_confidence;   //!< Confidence level
    int         m_max_iter;     //!< Maximum number of iterations
    bool        m_apply_edisp;  //!< Apply energy dispersion?
    bool        m_fit;          //!< Do fit?
    GEbounds    m_ebounds;      //!< Energy binning definition
    GFilename   m_outfile;      //!< Output ASCII file name
    GChatter    m_chatter;      //!< Chattiness

    // Protected members
    GMatrixSparse       m_covariance;      //!< Covariance matrix
    std::vector<double> m_energies;        //!< Energy values
    std::vector<double> m_intensities;     //!< Power law intensity
    std::vector<double> m_min_intensities; //!< Minimum intensities
    std::vector<double> m_max_intensities; //!< Maximum intensities
};


#endif /* CTBUTTERFLY_HPP */
