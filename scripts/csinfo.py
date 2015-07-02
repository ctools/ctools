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
import subprocess

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

   cscaldb       Lists available instrument response functions
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


def print_ctools_list():
    """Print available ctools with one-line description to the console.
    """
    print('\nAvailable ctools:')
    print(CTOOL_LIST)


def csinfo_setup_check():
    print('\nGammalib / ctools setup info:\n')

    all_ok = True

    log('   GAMMALIB environment variable ... ', end='')
    try:
        gammalib_env = os.environ['GAMMALIB']
        print('ok')
    except KeyError:
        print('NOT OK')
        all_ok = False

    log('   CTOOLS   environment variable ... ', end='')
    try:
        ctools_env = os.environ['CTOOLS']
        print('ok')
    except KeyError:
        print('NOT_OK')
        all_ok = False

    log('   gammalib Python import .......... ', end='')
    try:
        import gammalib
        print('ok')
    except:
        print('NOT OK')
        all_ok = False

    log('   ctools   Python import .......... ', end='')
    try:
        import ctools
        print('ok')
    except:
        print('NOT OK')
        all_ok = False

    if all_ok:
        print('\n   Your Gammalib / ctools setup is OK.\n')
    else:
        print('\n   WARNING: Your Gammalib / ctools setup is NOT OK!\n')
        print('      - Did you source gammalib-init.sh ?')
        print('      - Did you source gammalib-ctools.sh ?')
        print('      - Are you using the correct Python ?')

    print()

def csinfo_setup_info():

    log('   CTOOLS  environment variable ... ', end='')
    try:
        ctools_env = os.environ['CTOOLS']
        print('   $CTOOLS = {0}'.format(ctools_env))
    except KeyError:
        print('WARNING: $CTOOLS environment variable not set. Did you source ctools-init.sh?')
        all_ok = False

    try:
        import gammalib
        print('   Python `gammalib` path: {0}'.format(gammalib.__path__[0]))
    except:
        print('WARNING: Python `gammalib` import failed. Are you using the right Python?')

    try:
        import ctools
        print('   Python `ctools` path: {0}'.format(ctools.__path__[0]))
    except:
        print('WARNING: Python `ctools` import failed. Are you using the right Python?')

    print('   GAMMALIB CFLAGS: {0}'.format(get_pkg_config_info('cflags', 'gammalib')))
    print('   CTOOLS   CFLAGS: {0}'.format(get_pkg_config_info('cflags', 'ctools')))

    print()


def get_pkg_config_info(info, library):
    cmd = 'pkg-config --{0} {1}'.format(info, library)
    out = subprocess.getoutput(cmd)
    if 'was not found' in out:
        return 'Not available'
    else:
        return out


def csinfo_print_help():
    print('\nPrint info about Gammalib and ctools to the console.\n')
    print('Available commands:')
    print('   csinfo tools         List available ctools')
    print('   csinfo setup-check   Check Gammalib / ctools setup')
    print('   csinfo setup-info    Print Gammalib / ctools setup info')
    print()


def csinfo(argv):
    """Print info about Gammalib / ctools to the console.
    """
    if len(argv) <= 1:
        csinfo_print_help()
        sys.exit(0)

    cmd = argv[1]
    if cmd == 'list-tools':
        csinfo_list_tools()
    elif cmd == 'setup-check':
        csinfo_setup_check()
    elif cmd == 'setup-info':
        csinfo_setup_info()
    else:
        print('\nERROR: invalid command: `{0}`'.format(cmd))
        csinfo_print_help()
        sys.exit(-1)


if __name__ == '__main__':
    csinfo(sys.argv)
