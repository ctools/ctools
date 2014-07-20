/***************************************************************************
 *                     ctskymap - CTA sky mapping tool                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @file ctskymap.i
 * @brief CTA sky mapping tool definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "ctskymap.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class ctskymap
 *
 * @brief CTA sky mapping tool interface defintion
 ***************************************************************************/
class ctskymap : public GApplication  {
public:
    // Constructors and destructors
    ctskymap(void);
    explicit ctskymap(const GObservations& obs);
    ctskymap(int argc, char *argv[]);
    ctskymap(const ctskymap& app);
    virtual ~ctskymap(void);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void);
    GSkymap&       skymap(void);
    void           get_parameters(void);
    void           init_map(GCTAObservation* obs);
    void           map_events(GCTAObservation* obs);
};


/***********************************************************************//**
 * @brief CTA sky mapping tool Python extension
 ***************************************************************************/
%extend ctskymap {
    ctskymap copy() {
        return (*self);
    }
}
