/***************************************************************************
 *                     ctpsfcube - CTA PSF cube tool                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Chia-Chun Lu                                     *
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
 * @file ctpsfcube.i
 * @brief CTA PSF cube tool definition
 * @author Chia-Chun Lu
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctpsfcube.hpp"
#include "GTools.hpp"
%}

/***********************************************************************//**
 * @class ctpsfcube
 *
 * @brief CTA PSF cube tool
 ***************************************************************************/
class ctpsfcube : public GApplication  {

public:
    // Constructors and destructors
    ctpsfcube(void);
    explicit ctpsfcube(const GObservations& obs);
    ctpsfcube(int argc, char *argv[]);
    ctpsfcube(const ctpsfcube& app);
    virtual ~ctpsfcube(void);

    // Methods
    void                 clear(void);
    void                 execute(void);
    void                 run(void);
    void                 save(void);
    const GObservations& obs(void) const;
    const GCTAMeanPsf&   psfcube(void) const;
    void                 get_parameters(void);
    void                 get_obs(void);
    void                 set_response(void);
    void                 get_ebounds(void);
    void                 set_from_cntmap(const std::string& filename);
};

/***********************************************************************//**
 * @brief CTA PSF cube tool Python extension
 ***************************************************************************/
%extend ctpsfcube {
    ctpsfcube copy() {
        return (*self);
    }
}
