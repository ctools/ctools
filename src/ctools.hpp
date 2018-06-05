/***************************************************************************
 *    ctools.hpp - Cherenkov Telescope Array Analysis tools Header file    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file ctools.hpp
 * @brief ctools definitions
 * @author Juergen Knoedlseder
 */

#ifndef CTOOLS_HPP
#define CTOOLS_HPP


/***************************************************************************
 *                              Core services                              *
 ***************************************************************************/

/* __ Support functions __________________________________________________ */
#include "support.hpp"

/* __ Base classes _______________________________________________________ */
#include "ctool.hpp"
#include "ctobservation.hpp"
#include "ctlikelihood.hpp"


/***************************************************************************
 *                                   ctools                                *
 ***************************************************************************/

/* __ Simulation tools ___________________________________________________ */
#include "ctobssim.hpp"

/* __ Selection tools ____________________________________________________ */
#include "ctselect.hpp"

/* __ Binning tools ______________________________________________________ */
#include "ctbin.hpp"
#include "ctbkgcube.hpp"
#include "ctcubemask.hpp"
#include "ctedispcube.hpp"
#include "ctexpcube.hpp"
#include "ctpsfcube.hpp"

/* __ Imaging tools ______________________________________________________ */
#include "ctskymap.hpp"

/* __ Spectral tools _____________________________________________________ */

/* __ Timing tools _______________________________________________________ */
#include "ctphase.hpp"
#include "ctprob.hpp"

/* __ Likelihood tools ___________________________________________________ */
#include "ctbutterfly.hpp"
#include "cterror.hpp"
#include "ctlike.hpp"
#include "cttsmap.hpp"
#include "ctulimit.hpp"

/* __ Other tools ________________________________________________________ */
#include "ctmapcube.hpp"
#include "ctmodel.hpp"

#endif /* CTOOLS_HPP */
