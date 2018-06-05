/***************************************************************************
 *                      support - support functions                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file support
 * @brief Support function definitions
 * @author Juergen Knoedlseder
 */

#ifndef SUPPORT_HPP
#define SUPPORT_HPP

/* __ Includes ___________________________________________________________ */
#include "ctool.hpp"

/* __ Definitions ________________________________________________________ */

/* __ Prototypes _________________________________________________________ */
int  execute_ctool(ctool* tool);
void report_ctool_failure(const std::string& name, const std::string& message);

#endif /* SUPPORT_HPP */
