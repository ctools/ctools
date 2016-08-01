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
import sys
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


# ============================== #
# Get options from argument list #
# ============================== #
def get_arg_options(options, usage):
    """
    Get options from argument list

    Parameters
    ----------
    options : list of dict
        List of possible options and default values
    usage : str
        Usage string to be shown in case of a problem

    Returns
    -------
    options : list of dict
        Update list of options
    """
    # If there are no command line arguments then show usage string and
    # exit
    if len(sys.argv) < 1:
        print(usage)
        sys.exit()

    # First possible option is element 1
    i = 1

    # Loop over all arguments
    while i < len(sys.argv):

        # Search for options
        for option in options:
            if sys.argv[i] == option['option']:

                # If an option has been found then check if there is a
                # following parameter and extract that parameter as a
                # string
                if len(sys.argv) > i+1:
                    i += 1
                    try:
                        option['value'] = str(sys.argv[i])
                    except:
                        print(usage)
                        sys.exit()
                else:
                    print(usage)
                    sys.exit()

        # Next item
        i += 1

    # Return options
    return options
