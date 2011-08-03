/***************************************************************************
 *               ctobssim - CTA observation simulation tool                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file ctobssim.i
 * @brief CTA observation simulation tool Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctobssim.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class ctobssim
 *
 * @brief CTA data selection tool Python interface
 ***************************************************************************/
class ctobssim : public GApplication  {
public:
    // Constructors and destructors
    ctobssim(void);
    explicit ctobssim(GObservations obs);
    ctobssim(int argc, char *argv[]);
    ctobssim(const ctobssim& app);
    virtual ~ctobssim(void);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void);
    void           get_parameters(void);
    void           set_list(GCTAObservation* obs);
    void           simulate_source(GCTAObservation* obs, const GModels& models);
    void           simulate_background(GCTAObservation* obs, const GModels& models);
};


/***********************************************************************//**
 * @brief CTA observation simulation tool Python extension
 ***************************************************************************/
%extend ctobssim {
    ctobssim copy() {
        return (*self);
    }
}
