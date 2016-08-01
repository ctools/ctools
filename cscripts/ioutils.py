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


# ================ #
# Read pull values #
# ================ #
def read_pull_values(filename, parname):
    """
    Read pull values from pull distribution ASCII file
    
    Parameters
    ----------
    filename : str
        Pull distribution ASCII file
    parname : str
        Parameter

    Returns
    -------
    values : list of float
        List of pull values

    Raises
    ------
    NameError
    """
    # Initialise list of values
    values = []

    # Open reader
    reader = csv.reader(open(filename, 'r'), delimiter=',')

    # Read rows
    first = True
    index = -1
    for row in reader:

        # Get column index if first row
        if first:
            try:
                index = row.index(parname)
            except:
                print('ERROR: Parameter "'+parname+'" not found in list:')
                for p in row:
                    print('       "'+p+'"')
                raise NameError(parname)

        # Handle data rows
        else:
            values.append(float(row[index]))

        # Flag that first row has been passed
        first = False

    # Return values
    return values


# ============================================ #
# Get arguments and options from argument list #
# ============================================ #
def get_args_options(options, usage):
    """
    Get arguments and options from argument list

    Parameters
    ----------
    options : list of dict
        List of possible options and default values
    usage : str
        Usage string to be shown in case of a problem

    Returns
    -------
    args, options : tuple of list and list of dict
        Arguments and updated list of options
    """
    # Initialise list of arguments
    args = []

    # If there are no command line arguments then show usage string and
    # exit
    if len(sys.argv) < 1:
        print('Usage: %s' % usage)
        sys.exit()

    # First possible option is element 1
    i = 1

    # Loop over all arguments
    while i < len(sys.argv):

        # Initialise that no option was found
        option_found = False

        # Search for options
        for option in options:

            # Do we have an option
            if sys.argv[i] == option['option']:

                # If an option has been found then check if there is a
                # following parameter and extract that parameter as a
                # string. Quit in case that an exception occurs.
                if len(sys.argv) > i+1:
                    i += 1
                    try:
                        option['value'] = str(sys.argv[i])
                    except:
                        print('Usage: %s' % usage)
                        sys.exit()

                # ... there is no following parameter, hence write out usage
                # and quite
                else:
                    print('Usage: %s' % usage)
                    sys.exit()

                # We can break now as every option should only occur once
                option_found = True
                break

        # If no option was found then add the argument to list of arguments
        if not option_found:
            args.append(sys.argv[i])

        # Next item
        i += 1

    # Return arguments and options
    return args, options
