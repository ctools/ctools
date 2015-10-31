#!/usr/bin/env python
# ==========================================================================
# Print info about Gammalib / ctools to the console.
#
# Copyright (C) 2015 Christoph Deil
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


def get_command_output(cmd):
    """Utility function to get command output on old and new Python versions.
    """
    try:
        from subprocess import getoutput
    except ImportError:
        from commands import getoutput

    return getoutput(cmd)


# Copied and edited a bit from the `doc/source/reference_manual/index.rst` file
# The list there and here needs to be manually updated and kept in sync.
CTOOL_LIST = """
   ctbin         Generates counts cube
   ctbkgcube     Generates background cube
   ctbutterfly   Compute butterfly
   ctcubemask    Filter counts cube
   ctexpcube     Generates exposure cube
   ctlike        Performs maximum likelihood fitting
   ctmodel       Computes model counts cube
   ctobssim      Simulate CTA observations
   ctpsfcube     Generates point spread function cube
   ctselect      Selects event data
   ctskymap      Generates sky map
   cttsmap       Generates Test Statistics map
   ctulimit      Calculates upper limit
   cterror       Calculates likelihood profile errors
"""
CSCRIPT_LIST = """
   cscaldb       Lists available instrument response functions
   csinfo        Checks ctools and GammaLib installations
   csobsdef      Generates observation definition file
   cslightcrv    Computes lightcurve
   cspull        Generates pull distribution
   cssens        Computes CTA sensitivity
   csspec        Computes spectral points
   csresmap      Generates residual map
   cstsdist      Generates TS distribution
"""

def log(text, end='\n'):
    """Print text to console.

    Workaround to have an `end=''` option on old Pythons:
    http://stackoverflow.com/a/493399/498873
    """
    sys.stdout.write(text + end)


def csinfo_list_tools():
    """Print available ctools with one-line description to the console.
    """
    print('\nAvailable ctools:')
    print(CTOOL_LIST)
    print('Available cscripts:')
    print(CSCRIPT_LIST)


def csinfo_setup_check():
    print('\nGammalib / ctools setup check:\n')

    all_ok           = True
    gammalib_environ = False
    ctools_environ   = False
    gammalib_python  = False
    ctools_python    = False

    log('   GAMMALIB environment variable ... ', end='')
    try:
        gammalib_env = os.environ['GAMMALIB']
        print('ok')
    except KeyError:
        print('NOT OK')
        all_ok           = False
        gammalib_environ = True

    log('   CTOOLS   environment variable ... ', end='')
    try:
        ctools_env = os.environ['CTOOLS']
        print('ok')
    except KeyError:
        print('NOT_OK')
        all_ok         = False
        ctools_environ = True

    log('   gammalib Python import .......... ', end='')
    try:
        import gammalib
        print('ok')
    except:
        print('NOT OK')
        all_ok          = False
        gammalib_python = True

    log('   ctools   Python import .......... ', end='')
    try:
        import ctools
        print('ok')
    except:
        print('NOT OK')
        all_ok        = False
        ctools_python = True

    if all_ok:
        print('\n   ===> Your Gammalib / ctools setup is OK.')
    else:
        print('\n   ===> WARNING: Your Gammalib / ctools setup is NOT OK!\n')
        if gammalib_environ:
            print('      - Have you set the GAMMALIB environment variable?')
        if ctools_environ:
            print('      - Have you set the CTOOLS environment variable?')
        if gammalib_python or ctools_python:
            print('      - Did you source gammalib-init.sh?')
            print('      - Did you source ctools-init.sh ?')
            print('      - Are you using the correct Python ?')

        print('')
        print('   Gammalib:  http://cta.irap.omp.eu/gammalib-devel/')
        print('   ctools:    http://cta.irap.omp.eu/ctools-devel/')

    print('')


def csinfo_setup_info():
    print('\nGammalib / ctools setup info:\n')

    gammalib_version = 'TODO'  # a function needs to be added to the library to get the version
    print('   Gammalib version ................. {0}'.format(gammalib_version))

    ctools_version = 'TODO'  # a function needs to be added to the library to get the version
    print('   ctools   version ................. {0}'.format(ctools_version))

    gammalib_env = os.environ.get('GAMMALIB', 'Not found')
    print('   $GAMMALIB environment variable ... {0}'.format(gammalib_env))

    ctools_env = os.environ.get('CTOOLS', 'Not found')
    print('   $CTOOLS   environment variable ... {0}'.format(ctools_env))

    print('   Python executable ................ {0}'.format(sys.executable))

    log('   gammalib  Python package ......... ', end='')
    try:
        import gammalib
        print(gammalib.__path__[0])
    except ImportError:
        print('Not found')

    log('   ctools    Python package  ........ ', end='')
    try:
        import ctools
        print(ctools.__path__[0])
    except ImportError:
        print('Not found')

    print('   GAMMALIB  CFLAGS ................. {0}'.format(get_pkg_config_info('cflags', 'gammalib')))
    print('   CTOOLS    CFLAGS ................. {0}'.format(get_pkg_config_info('cflags', 'ctools')))

    print('   GAMMALIB  LIBS   ................. {0}'.format(get_pkg_config_info('libs', 'gammalib')))
    print('   CTOOLS    LIBS   ................. {0}'.format(get_pkg_config_info('libs', 'ctools')))

    print('')


def get_pkg_config_info(info, library):
    cmd = 'pkg-config --{0} {1}'.format(info, library)
    out = get_command_output(cmd)
    if 'was not found' in out:
        return 'Not available'
    else:
        return out


def csinfo_print_help():
    print('\nPrint info about Gammalib and ctools to the console.\n')
    print('Available commands:')
    print('   csinfo list    List available ctools and cscripts')
    print('   csinfo check   Check Gammalib / ctools setup')
    print('   csinfo info    Print Gammalib / ctools setup info')
    print('')


def csinfo(argv):
    """Print info about Gammalib / ctools to the console.
    """
    if len(argv) <= 1:
        csinfo_print_help()
        sys.exit(0)

    cmd = argv[1]
    if cmd == 'list':
        csinfo_list_tools()
    elif cmd == 'check':
        csinfo_setup_check()
    elif cmd == 'info':
        csinfo_setup_info()
    else:
        print('\nERROR: invalid command: `{0}`'.format(cmd))
        csinfo_print_help()
        sys.exit(-1)


if __name__ == '__main__':
    csinfo(sys.argv)
