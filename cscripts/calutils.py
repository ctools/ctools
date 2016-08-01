# ==========================================================================
# Utility functions for calibration handling
#
# Copyright (C) 2016 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import gammalib


# ======================= #
# Create CIF binary table #
# ======================= #
def create_cif_table():
    """
    Create Calibration Database Index File binary table
    """
    # Create binary table
    table = gammalib.GFitsBinTable()

    # Append columns. Reference: CAL/GEN/92-008
    table.append(gammalib.GFitsTableStringCol("TELESCOP", 0, 10))
    table.append(gammalib.GFitsTableStringCol("INSTRUME", 0, 10))
    table.append(gammalib.GFitsTableStringCol("DETNAM", 0, 20))
    table.append(gammalib.GFitsTableStringCol("FILTER", 0, 10))
    table.append(gammalib.GFitsTableStringCol("CAL_DEV", 0, 20))
    table.append(gammalib.GFitsTableStringCol("CAL_DIR", 0, 70))
    table.append(gammalib.GFitsTableStringCol("CAL_FILE", 0, 40))
    table.append(gammalib.GFitsTableStringCol("CAL_CLAS", 0, 3))
    table.append(gammalib.GFitsTableStringCol("CAL_DTYP", 0, 4))
    table.append(gammalib.GFitsTableStringCol("CAL_CNAM", 0, 20))
    table.append(gammalib.GFitsTableStringCol("CAL_CBD", 0, 70, 9))
    table.append(gammalib.GFitsTableShortCol("CAL_XNO", 0))
    table.append(gammalib.GFitsTableStringCol("CAL_VSD", 0, 10))
    table.append(gammalib.GFitsTableStringCol("CAL_VST", 0, 8))
    table.append(gammalib.GFitsTableDoubleCol("REF_TIME", 0))
    table.append(gammalib.GFitsTableShortCol("CAL_QUAL", 0))
    table.append(gammalib.GFitsTableStringCol("CAL_DATE", 0, 8))
    table.append(gammalib.GFitsTableStringCol("CAL_DESC", 0, 70))

    # Set keywords. Reference: CAL/GEN/92-008
    table.extname("CIF")
    table.card("CIFVERSN", "1992a", "Version of CIF format")

    # Return table
    return table
