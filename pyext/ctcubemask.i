/***************************************************************************
 *                    ctcubemask - CTA mask cube tool                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 Chia-Chun Lu                                        *
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
 * @file ctcubemask.i
 * @brief CTA data selection tool Python interface definition
 * @author Chia-Chun Lu
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctcubemask.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class ctcubemask
 *
 * @brief CTA data selection tool Python interface
 ***************************************************************************/
class ctcubemask : public GApplication  {
public:
    // Constructors and destructors
    ctcubemask(void);
    explicit ctcubemask(GObservations obs);
    ctcubemask(int argc, char *argv[]);
    ctcubemask(const ctcubemask& app);
    virtual ~ctcubemask(void);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           apply_mask(GCTAObservation* obs, const std::string& filename);
};


/***********************************************************************//**
 * @brief CTA data selection tool Python extension
 ***************************************************************************/
%extend ctcubemask {
    ctcubemask copy() {
        return (*self);
    }
}
