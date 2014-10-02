/***************************************************************************
 *                      cttsmap - CTA data binning tool                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file cttsmap.i
 * @brief CTA data binning tool Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "cttsmap.hpp"
#include "GTools.hpp"
%}


/***********************************************************************
//**
* @class cttsmap
*
* @brief CTA data binning tool Python interface
***************************************************************************/
class cttsmap : public GApplication  {
public:
    // Constructors and destructors
    cttsmap(void);
    explicit cttsmap(const GObservations& obs);
    cttsmap(int argc, char *argv[]);
    cttsmap(const cttsmap& app);
    virtual ~cttsmap(void);
    
    // Methods
    void                 clear(void);
    void                 execute(void);
    void                 run(void);
    void                 save(void);
    const GObservations& obs(void) const;
    const GSkymap& tsmap(void) const;
    const GSkymap& fluxmap(void) const;
    const GSkymap& indexmap(void) const;
    const GSkymap& statusmap(void) const;

    
};


/***********************************************************************
 //**
* @brief CTA data binning tool Python extension
***************************************************************************/
%extend cttsmap {
    cttsmap copy() {
        return (*self);
    }
}
