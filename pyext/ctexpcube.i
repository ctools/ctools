/***************************************************************************
 *                   ctexpcube - CTA exposure cube tool                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file ctexpcube.i
 * @brief CTA exposure cube tool definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctexpcube.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class ctexpcube
 *
 * @brief CTA exposure cube tool
 ***************************************************************************/
class ctexpcube : public GApplication  {

public:
    // Constructors and destructors
    ctexpcube(void);
    explicit ctexpcube(GObservations obs);
    ctexpcube(int argc, char *argv[]);
    ctexpcube(const ctexpcube& app);
    virtual ~ctexpcube(void);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void);
    void           get_parameters(void);
    void           get_obs(void);
    void           set_response(void);
    void           get_ebounds(void);
    void           set_from_cntmap(const std::string& filename);
};


/***********************************************************************//**
 * @brief CTA exposure map tool Python extension
 ***************************************************************************/
%extend ctexpcube {
    ctexpcube copy() {
        return (*self);
    }
}
