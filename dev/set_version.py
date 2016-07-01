#! /usr/bin/env python
# ==========================================================================
# Set the version number of the software
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
import sys


# ========================== #
# Get current version number #
# ========================== #
def get_current():
    """
    Get current version number

    Returns
    -------
    version : str
        Version number
    """
    # Initialse version number of empty string
    version = ''

    # Open configure.ac file
    f = open('configure.ac', 'r')

    # Read all lines
    for line in f:

        # Search for pattern
        if line.startswith('AC_INIT('):

            # Extract version number
            start   = line.find(',')+1
            stop    = line.find(',', start)
            version = line[start:stop].strip().lstrip('[').rstrip(']')

    # Close file
    f.close()

    # Return version
    return version


# ============================== #
# Set version number in one file #
# ============================== #
def set_version(filename, version, current):
    """
    Set version number

    Parameters
    ----------
    filename : str
        File name
    version : str
        Version number
    current : str
        Current version number
    """
    # Open file in read mode
    f = open(filename, 'r')

    # Read file content
    content = f.read()

    # Close file
    f.close()

    # Replace version number
    new_content = content.replace(current, version)
    #print(new_content)

    # Open file in write mode
    f = open(filename, 'w')

    # Write file content
    f.write(new_content)

    # Close file
    f.close()

    # Return
    return


# =============================== #
# Set version number in all files #
# =============================== #
def set_versions(version, current):
    """
    Set version numbers in all files

    Parameters
    ----------
    version : str
        Version number
    current : str
        Current version number
    """
    # List of all relevant files, including their access paths from the
    # source directory
    filenames = ['configure.ac',
                 'README.md',
                 'sonar-project.properties',
                 'doc/Doxyfile',
                 'doc/source/conf.py']

    # Set versions by loop over all files
    for filename in filenames:
        set_version(filename, version, current)

    # Return
    return


# =============================================== #
# Set special version number in configure.ac file #
# =============================================== #
def set_special_configure(version, current):
    """
    Set special version numbers in configure.ac file

    Parameters
    ----------
    version : str
        Version number
    current : str
        Current version number
    """
    # Generate current version number
    split = current.split('.')
    cur   = '%s:%s:%s' % (split[0], split[1], split[2])

    # Generate new version number
    split = version.split('.')
    ver   = '%s:%s:%s' % (split[0], split[1], split[2])

    # Set version number
    set_version('configure.ac', ver, cur)

    # Return
    return


# ===================================================== #
# Set special version number in doc/source/conf.py file #
# ===================================================== #
def set_special_conf(version, current):
    """
    Set special version numbers in doc/source/conf.py file

    Parameters
    ----------
    version : str
        Version number
    current : str
        Current version number
    """
    # Generate current version number
    split = current.split('.')
    cur   = '%s.%s' % (split[0], split[1])

    # Generate new version number
    split = version.split('.')
    ver   = '%s.%s' % (split[0], split[1])

    # Set version number
    set_version('doc/source/conf.py', ver, cur)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Check for version number
    if len(sys.argv) < 2:
        sys.exit('Usage: set_version.py version')

    # Extract version number
    version = sys.argv[1]

    # Get the current version
    current = get_current()

    # Set versions
    set_versions(version, current)

    # Set special version number:
    set_special_configure(version, current) # CTOOLS_LT_VERSION="x:y:z"
    set_special_conf(version, current)      # version = 'x.y'

