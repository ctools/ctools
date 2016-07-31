# ==========================================================================
# Utility functions for input and output
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
import csv
import gammalib
import ctools


# ========================= #
# Write one row of CSV file #
# ========================= #
def write_csv_row(outfile, row, colnames, colvalues):
    """
    Write one row into CSV file

    Parameters
    ----------
    outfile : str
        Output file name
    row : int
        Row number
    colnames : list of str
        Column names
    colvalues : list of float
        Column values
    """
    # If this is the first row then open a new file and write the header
    if row == 0:
        f       = open(outfile, 'w')
        writer  = csv.DictWriter(f, colnames)
        headers = {}
        for colname in colnames:
            headers[colname] = colname
        writer.writerow(headers)

    # ... otherwise append to an existing file
    else:
        f = open(outfile, 'a')

    # Write out row
    writer = csv.DictWriter(f, colnames)
    writer.writerow(colvalues)

    # Close file
    f.close()

    # Return
    return
