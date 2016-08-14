#!/usr/bin/env python
# ==========================================================================
# Print info about Gammalib / ctools to the console
#
# Copyright (C) 2015-2016 Christoph Deil
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
import os

# List of available ctools. The list was copied from
# `doc/source/users/reference_manual/index.rst`
# and needs to be manually updated and kept in sync
CTOOL_LIST = """
   ctbin         Generates counts cube
   ctbkgcube     Generates background cube
   ctbutterfly   Compute butterfly
   ctcubemask    Filter counts cube
   ctedispcube   Generates energy dispersion cube
   cterror       Calculates likelihood profile errors
   ctexpcube     Generates exposure cube
   ctlike        Performs maximum likelihood fitting
   ctmapcube     Generates a map cube
   ctmodel       Computes model counts cube
   ctobssim      Simulate CTA observations
   ctpsfcube     Generates point spread function cube
   ctselect      Selects event data
   ctskymap      Generates sky map
   cttsmap       Generates Test Statistics map
   ctulimit      Calculates upper limit
"""

# List of available cscripts. The list was copied from
# `doc/source/users/reference_manual/index.rst`
# and needs to be manually updated and kept in sync
CSCRIPT_LIST = """
   cscaldb       Lists available instrument response functions
   cslightcrv    Computes lightcurve
   csinfo        Checks ctools and GammaLib installations
   csmodelinfo   Shows model container content
   csmodelmerge  Merges several model containers into one file
   csobsdef      Generates observation definition file
   csobsinfo     Shows observation container content
   cspull        Generates pull distribution
   csresmap      Generates residual map
   cssens        Computes CTA sensitivity
   csspec        Computes spectral points
   cstsdist      Generates TS distribution
   cstsmapsplit  Creates commands to split the Test Statistic map computations
   cstsmapmerge  Merges slices from Test Statistic map computations\n
   csobs2caldb   Creates a caldb entry from an input observation
   csroot2caldb  Creates a caldb entry from a ROOT performance file
   csiactdata    Shows information about IACT data available on the user machine
   csiactobs     Generates observation definition file for IACT data from observation IDs
   csfindobs     Generates a list of IACT observation IDs
   csiactcopy    Copies IACT data from one location to another
"""


# ============================================================ #
# Workaround for command output on old and new Python versions #
# ============================================================ #
def get_command_output(cmd):
    """
    Utility function to get command output on old and new Python versions

    Parameters
    ----------
    cmd : str
        Command string

    Returns
    -------
    output : str
        Command output
    """
    # Try to import getoutput from the "subprocess" module, alternatively get
    # it from the "commands" module
    try:
        from subprocess import getoutput
    except ImportError:
        from commands import getoutput

    # Return command output
    return getoutput(cmd)


# ====================== #
# Print list information #
# ====================== #
def csinfo_list_tools():
    """
    Print available ctools with one-line description to the console
    """
    # Print list of available ctools
    print('\nAvailable ctools:')
    print(CTOOL_LIST)

    # Print list of available cscripts
    print('Available cscripts:')
    print(CSCRIPT_LIST)

    # Return
    return


# ======================= #
# Print check information #
# ======================= #
def csinfo_setup_check():
    """
    Print check information
    """
    # Print header
    print('\nGammalib / ctools setup check:\n')

    # Check if the "GAMMALIB" environment variable is set
    sys.stdout.write('   GAMMALIB environment variable ... ')
    gammalib_environ = 'GAMMALIB' in os.environ
    if gammalib_environ:
        print('ok')
    else:
        print('NOT OK')

    # Check if the "CTOOLS" environment variable is set
    sys.stdout.write('   CTOOLS   environment variable ... ')
    ctools_environ = 'CTOOLS' in os.environ
    if ctools_environ:
        print('ok')
    else:
        print('NOT OK')

    # Check if gammalib Python module import works
    sys.stdout.write('   gammalib Python import .......... ')
    try:
        import gammalib
        print('ok')
        gammalib_python = True
    except:
        print('NOT OK')
        gammalib_python = False

    # Check if ctools Python module import works
    sys.stdout.write('   ctools   Python import .......... ')
    try:
        import ctools
        print('ok')
        ctools_python = True
    except:
        print('NOT OK')
        ctools_python = False

    # Check if cscripts Python module import works
    sys.stdout.write('   cscripts Python import .......... ')
    try:
        import cscripts
        print('ok')
        cscripts_python = True
    except:
        print('NOT OK')
        cscripts_python = False

    # Set flag that everything is okay
    all_ok = gammalib_environ and ctools_environ and gammalib_python and \
             ctools_python and cscripts_python

    # Signal if everything is okay, or notify about the problems that have
    # been encountered
    if all_ok:
        print('\n   ===> Your Gammalib / ctools setup is OK.')
    else:
        print('\n   ===> WARNING: Your Gammalib / ctools setup is NOT OK!\n')
        if not gammalib_environ:
            print('      - Have you set the GAMMALIB environment variable?')
        if not ctools_environ:
            print('      - Have you set the CTOOLS environment variable?')
        if not gammalib_python or not ctools_python or not cscripts_python:
            print('      - Did you source gammalib-init.sh?')
            print('      - Did you source ctools-init.sh ?')
            print('      - Are you using the correct Python ?')
        print('')
        print('   Gammalib:  http://cta.irap.omp.eu/gammalib/')
        print('   ctools:    http://cta.irap.omp.eu/ctools/')
    print('')

    # Return
    return


# ======================= #
# Print setup information #
# ======================= #
def csinfo_setup_info():
    """
    Print setup information
    """
    # Try importing GammaLib, and flag if this was unsuccessful
    try:
        import gammalib
        gammalib_path    = gammalib.__path__[0]
        gammalib_version = gammalib.__version__
    except ImportError:
        gammalib_path    = 'Not found'
        gammalib_version = 'Unknown'

    # Try importing ctools, and flag if this was unsuccessful
    try:
        import ctools
        ctools_path    = ctools.__path__[0]
        ctools_version = ctools.__version__
    except ImportError:
        ctools_path    = 'Not found'
        ctools_version = 'Unknown'

    # Try importing cscripts, and flag if this was unsuccessful
    try:
        import cscripts
        cscripts_path    = cscripts.__path__[0]
        cscripts_version = cscripts.__version__
    except ImportError:
        cscripts_path    = 'Not found'
        cscripts_version = 'Unknown'

    # Get GAMMALIB and CTOOLS environment variables
    gammalib_env = os.environ.get('GAMMALIB', 'Not found')
    ctools_env   = os.environ.get('CTOOLS', 'Not found')

    # Print setup information
    print('\nGammalib / ctools setup info:\n')
    print('   Gammalib  version ................ %s' % gammalib_version)
    print('   ctools    version ................ %s' % ctools_version)
    print('   cscripts  version ................ %s' % cscripts_version)
    print('   $GAMMALIB environment variable ... %s' % gammalib_env)
    print('   $CTOOLS   environment variable ... %s' % ctools_env)
    print('   Python executable ................ %s' % sys.executable)
    print('   gammalib  Python module .......... %s' % gammalib_path)
    print('   ctools    Python module .......... %s' % ctools_path)
    print('   cscripts  Python module .......... %s' % cscripts_path)
    print('   GAMMALIB  CFLAGS ................. %s' % get_pkg_config_info('cflags', 'gammalib'))
    print('   CTOOLS    CFLAGS ................. %s' % get_pkg_config_info('cflags', 'ctools'))
    print('   GAMMALIB  LIBS   ................. %s' % get_pkg_config_info('libs', 'gammalib'))
    print('   CTOOLS    LIBS   ................. %s' % get_pkg_config_info('libs', 'ctools'))
    print('')

    # Return
    return


# =============================== #
# Get information from pkg-config #
# =============================== #
def get_pkg_config_info(info, library):
    """
    Get information from pkg-config

    Parameters
    ----------
    info : str
        Information to extract
    library : str
        Library name for which information is to be extracted

    Returns
    -------
    out : str
        Information string, 'Not available' if pkg-config is not installed
    """
    # Set pkg-config command
    cmd = 'pkg-config --%s %s' % (info, library)

    # Execute command
    out = get_command_output(cmd)

    # If the command returns "not found" then pkg-config is not available
    if 'not found' in out:
        out = 'Not available'

    # Return output
    return out


# ============================ #
# Print help text into console #
# ============================ #
def csinfo_print_help():
    """
    Print help text into console
    """
    # Print help text
    print('\nPrint info about Gammalib and ctools to the console.\n')
    print('Available commands:')
    print('   csinfo list    List available ctools and cscripts')
    print('   csinfo check   Check Gammalib / ctools setup')
    print('   csinfo info    Print Gammalib / ctools setup info')
    print('')

    # Return
    return


# ============================ #
# Print information to console #
# ============================ #
def csinfo(argv):
    """
    Print information about Gammalib and ctools to the console
    """
    # If there are no command line arguments then print help text and exit
    # with success
    if len(argv) <= 1:
        csinfo_print_help()
        sys.exit(0)

    # Dispatch according to the command line argument
    cmd = argv[1]
    if cmd == 'list':
        csinfo_list_tools()
    elif cmd == 'check':
        csinfo_setup_check()
    elif cmd == 'info':
        csinfo_setup_info()
    else:
        print('\nERROR: invalid command: "%s"' % cmd)
        csinfo_print_help()
        sys.exit(-1)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Run the script
    csinfo(sys.argv)
