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
 * @brief Support function implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <cstdlib>         // std::getenv() function
//#include <cstdio>          // std::fopen(), etc. functions
#include "GammaLib.hpp"
#include "GCTALib.hpp"
#include "support.hpp"
//#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Coding definitions _________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                            Support functions                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Execute ctool
 *
 * @param[in] tool Pointer to ctool.
 * @return Return code
 *
 * Executes a ctool. The execution is enclosed in a try-catch statement so
 * that any exceptions that occur during execution will be properly logged
 * and written into the standard output.
 ***************************************************************************/
int execute_ctool(ctool* tool)
{
    // Initialise return code
    int rc = 1;

    // Try executing tool. This catches any exceptions that occur during
    // the execution.
    try {

        // Execute ctool
        tool->execute();

        // Signal success
        rc = 0;

    } // endtry

    // Handle any exceptions that have occured
    catch (std::exception &e) {

        // Extract error message
        std::string message = e.what();
        std::string signal  = "*** ERROR encounterted in the execution of "+
                              tool->name()+". Run aborted ...";

        // Write error in logger
        tool->log << signal  << std::endl;
        tool->log << message << std::endl;

        // Write error in standard output
        std::cout << signal  << std::endl;
        std::cout << message << std::endl;

    } // endcatch: catched any exception

    // Return return code
    return rc;
}


/***********************************************************************//**
 * @brief Report ctool failure
 *
 * @param[in] name Name of ctool.
 * @param[in] message Exception message.
 *
 * Reports a ctools failure.
 ***************************************************************************/
void report_ctool_failure(const std::string& name, const std::string& message)
{
    // Set signal string
    std::string signal  = "*** ERROR encounterted in the execution of "+
                          name+". Run aborted ...";

    // Write error in standard output
    std::cout << signal  << std::endl;
    std::cout << message << std::endl;

    // Return
    return;
}
