/***************************************************************************
 *                     ctmodel - CTA counts model tool                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file ctmodel.i
 * @brief CTA counts model tool definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctmodel.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class ctmodel
 *
 * @brief CTA counts model tool interface defintion
 ***************************************************************************/
class ctmodel : public GApplication  {

public:
    // Constructors and destructors
    ctmodel(void);
    explicit ctmodel(GObservations obs);
    ctmodel(int argc, char *argv[]);
    ctmodel(const ctmodel& app);
    virtual ~ctmodel(void);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           setup_obs(void);
    void           model_map(GCTAObservation* obs, const GModels& models);
};


/***********************************************************************//**
 * @brief CTA counts model tool Python extension
 ***************************************************************************/
%extend ctmodel {
    ctmodel copy() {
        return (*self);
    }
}
