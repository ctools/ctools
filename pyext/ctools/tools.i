/***************************************************************************
 *                           ctools - SWIG file                            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * ----------------------------------------------------------------------- *
 * Usage:                                                                  *
 * swig -c++ -python -Wall tools.i                                         *
 ***************************************************************************/
/**
 * @file tools.i
 * @brief ctools SWIG file
 * @author Juergen Knoedlseder
 */
%module tools
%feature("autodoc", "1");

/* __ Headers needed for compilation _____________________________________ */
%{
#include <stddef.h>
%}

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i
%include "std_vector.i"
namespace std {
   %template(vector_int)    vector<int>;
   %template(vector_double) vector<double>;
   %template(vector_bool)   vector<bool>;
};

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i";

/* __ Inform about base classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";
%import(module="gammalib.app")  "GApplication.i";

/* __ ctools base classes ________________________________________________ */
%include "ctool.i"
%include "ctobservation.i"
%include "ctlikelihood.i"

/* __ ctools _____________________________________________________________ */
%include "ctbin.i"
%include "ctobssim.i"
%include "ctlike.i"
%include "ctmodel.i"
%include "ctselect.i"
%include "ctskymap.i"
%include "ctexpcube.i"
%include "ctpsfcube.i"
%include "ctedispcube.i"
%include "ctbkgcube.i"
%include "ctmapcube.i"
%include "ctcubemask.i"
%include "cttsmap.i"
%include "ctbutterfly.i"
%include "ctulimit.i"
%include "cterror.i"

/* __ Test function ______________________________________________________ */
%pythoncode %{
def test():
    from ctools.tests import test_python_ctools
    test_python_ctools.test(installed=True)
%}
