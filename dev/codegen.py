#! /usr/bin/env python
# ==========================================================================
# ctools code generator
#
# Copyright (C) 2018 Juergen Knoedlseder
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
from __future__ import print_function
import os
import sys
import fileinput
from datetime import date
from builtins import input
import gammalib


# ================= #
# Confirm something #
# ================= #
def confirm(text):
    """
    Confirm something

    Parameters
    ----------
    text : str
        Thing to confirm

    Returns
    -------
    flag : bool
        True if something is confirmed
    """
    # Confirmation loop
    waiting = True
    while waiting:
        confirmation = str(input(text+' (y/n): '))
        if confirmation == 'q':
            sys.exit()
        elif confirmation == 'y':
            flag    = True
            waiting = False
        elif confirmation == 'n':
            flag    = False
            waiting = False

    # Return confirmation flag
    return flag


# ============ #
# Get response #
# ============ #
def response(text, confirm_response=False):
    """
    Get response

    Parameters
    ----------
    text : str
        Thing we want to get a response to
    confirm_response : bool
        True if response should be confirmed

    Returns
    -------
    response : str
        Reponse
    """
    # Looping until response is accepted
    looping = True
    while looping:

        # Get response from input
        waiting = True
        while waiting:
            response = str(input('%s: ' % text))
            if response == 'q':
                sys.exit()
            else:
                waiting = False

        # If the version number change is confirmed then do the change
        if confirm_response:
            if confirm('Is "%s" ok?' % response):
                looping = False
        else:
            looping = False

    # Return response
    return response


# ======== #
# Add file #
# ======== #
def add_file(infile, outfile, tokens):
    """
    Add a file, replacing all tokens

    Parameters
    ----------
    infile : str
        Input filename
    outfile : str
        Output filename
    tokens : list of dict
        Token list
    """
    # Create output file
    file = open(outfile, "w")

    # Initialise number of header lines
    n_header = 0

    # Loop over input file
    for line in open(infile, 'r'):

        # Signal if this is a header line
        is_header_ast  = False
        is_header_hash = False
        if len(line) > 75:
            is_header_ast  = line[1] == '*' and line[75] == '*'
            is_header_hash = line[0] == '#' and line[1] == ' ' and line[76] == '#'
            n_header += 1

        # Replace tokens in line
        for token in tokens:
            line = line.replace(token['pattern'], token['string'])

        # If we have a header line than format the line
        if is_header_ast:
            text = line.strip(' *\n')+'                            '
            if n_header == 2:
                line = ' *'+gammalib.centre(text,73)+'*\n'
            elif n_header == 4:
                line = ' *'+gammalib.left('  '+text,73)+'*\n'
        if is_header_hash:
            text = line.strip('#\n')+'                            '
            line = '#'+text[0:75]+'#\n'

        # Write out line
        file.write(line)

    # Close file
    file.close()

    # Return
    return


# ================ #
# Set ctool tokens #
# ================ #
def set_ctool_tokens(classname, basename, what, author, email, affiliation):
    """
    Set replacement tokens for a ctool
    """
    # Get current year
    year = str(date.today().year)

    # Set tokens
    tokens = [{'pattern': 'ctool_%s' % basename, 'string': classname},
              {'pattern': 'CTOOL_%s' % basename.upper(), 'string': classname.upper()},
              {'pattern': 'xxx', 'string': classname},
              {'pattern': 'XXX', 'string': classname.upper()},
              {'pattern': '[WHAT]', 'string': what},
              {'pattern': '[what]', 'string': what.lower()},
              {'pattern': '[AUTHOR]', 'string': author},
              {'pattern': '[EMAIL]', 'string': email},
              {'pattern': '[AFFILIATION]', 'string': affiliation},
              {'pattern': '[YEAR]', 'string': year}]

    # Return tokens
    return tokens


# ================== #
# Set cscript tokens #
# ================== #
def set_cscript_tokens(classname, basename, what, author, email, affiliation):
    """
    Set replacement tokens for a cscript
    """
    # Get current year
    year = str(date.today().year)

    # Set tokens
    tokens = [{'pattern': 'cscript_%s' % basename, 'string': classname},
              {'pattern': 'CSCRIPT_%s' % basename.upper(), 'string': classname.upper()},
              {'pattern': 'xxx', 'string': classname},
              {'pattern': 'XXX', 'string': classname.upper()},
              {'pattern': '[WHAT]', 'string': what},
              {'pattern': '[what]', 'string': what.lower()},
              {'pattern': '[AUTHOR]', 'string': author},
              {'pattern': '[EMAIL]', 'string': email},
              {'pattern': '[AFFILIATION]', 'string': affiliation},
              {'pattern': '[YEAR]', 'string': year}]

    # Return tokens
    return tokens


# ========= #
# Add ctool #
# ========= #
def add_ctool(name, tokens, baseclass):
    """
    Add a ctool

    Parameters
    ----------
    name : str
        ctool name
    tokens : str
        Tokens for replacement
    baseclass : str
        Baseclass
    """
    # Create directory structure
    if not os.path.isdir('src/%s' % name):
        os.mkdir('src/%s' % name)

    # Set template file names
    inctemp  = 'src/template/ctool_%s.hpp' % (baseclass)
    srctemp  = 'src/template/ctool_%s.cpp' % (baseclass)
    pytemp   = 'src/template/ctool_%s.i'   % (baseclass)
    maintemp = 'src/template/main_%s.cpp'  % (baseclass)
    partemp  = 'src/template/ctool.par'
    maketemp = 'src/template/ctool_Makefile.am'
    testtemp = 'src/template/ctool_test.py'
    doctemp  = 'src/template/ctool.rst'

    # Set destination file names
    incfile  = 'src/%s/%s.hpp'      % (name, name)
    srcfile  = 'src/%s/%s.cpp'      % (name, name)
    pyfile   = 'pyext/%s.i'         % (name)
    mainfile = 'src/%s/main.cpp'    % (name)
    parfile  = 'src/%s/%s.par'      % (name, name)
    makefile = 'src/%s/Makefile.am' % (name)
    testfile = 'test/test_%s.py'    % (name)
    docfile  = 'doc/source/users/reference_manual/%s.rst' % (name)
    
    # Add files
    add_file(inctemp,  incfile,  tokens)
    add_file(srctemp,  srcfile,  tokens)
    add_file(partemp,  parfile,  tokens)
    add_file(maintemp, mainfile, tokens)
    add_file(maketemp, makefile, tokens)
    add_file(pytemp,   pyfile,   tokens)
    add_file(testtemp, testfile, tokens)
    add_file(doctemp,  docfile,  tokens)

    # Update src/Makefile.am
    filename    = 'src/Makefile.am'
    insertline1 = '          %s ' % (name)
    insertline2 = '                      %s/lib%s.la ' % (name, name)
    for line in fileinput.FileInput(filename,inplace=1):
        if 'SUBDIRS' in line.replace(' ', ''):
            last = line[-2]
            line = line.replace(line, line+insertline1+last+'\n')
        if 'libctools_la_LIBADD' in line.replace(' ', ''):
            last = line[-2]
            line = line.replace(line, line+insertline2+last+'\n')
        print(line,end=''),

    # Update Python module
    for line in fileinput.FileInput('pyext/ctools/tools.i',inplace=1):
        if '%include "ctbin.i"' in line:
            print('%%include "%s.i"' % (name))
        print(line,end=''),

    # Update unit test
    for line in fileinput.FileInput('test/test_python_ctools.py',inplace=1):
        if 'import test_ctbin' in line:
            print('import test_%s' % (name))
        elif 'test_ctbin.Test()' in line:
            print('             test_%s.Test(),' % (name))
        print(line,end=''),

    # Update export of unit test
    for line in fileinput.FileInput('pyext/Makefile.am',inplace=1):
        if '$(top_srcdir)/test/test_ctbin.py' in line:
            print('              $(top_srcdir)/test/test_%s.py \\' % (name))
        print(line,end=''),

    # Update configure.ac
    for line in fileinput.FileInput('configure.ac',inplace=1):
        if 'src/ctbin/Makefile' in line:
            print('                 src/%s/Makefile' % (name))
        print(line,end=''),

    # Update reference manual index file
    for line in fileinput.FileInput('doc/source/users/reference_manual/index.rst',inplace=1):
        if 'ctbin --- Generates counts cube <ctbin>' in line:
            print('   %s --- ToDo: Describe what tool is doing <%s>' % (name, name))
        print(line,end=''),

    # Return
    return


# =========== #
# Add cscript #
# =========== #
def add_cscript(name, tokens, baseclass):
    """
    Add a cscript

    Parameters
    ----------
    name : str
        cscript name
    tokens : str
        Tokens for replacement
    baseclass : str
        Baseclass
    """
    # Set template file names
    pytemp   = 'src/template/cscript_%s.py' % (baseclass)
    partemp  = 'src/template/cscript.par'
    testtemp = 'src/template/cscript_test.py'
    doctemp  = 'src/template/cscript.rst'

    # Set destination file names
    pyfile   = 'cscripts/%s.py'  % (name)
    parfile  = 'cscripts/%s.par' % (name)
    testfile = 'test/test_%s.py' % (name)
    docfile  = 'doc/source/users/reference_manual/%s.rst' % (name)
    
    # Add files
    add_file(pytemp,   pyfile,   tokens)
    add_file(partemp,  parfile,  tokens)
    add_file(testtemp, testfile, tokens)
    add_file(doctemp,  docfile,  tokens)

    # Set file permissions
    os.chmod(pyfile, 0o755)

    # Update cscripts/__init__.py.in file
    for line in fileinput.FileInput('cscripts/__init__.py.in',inplace=1):
        if '"csworkflow",' in line:
            print('    "%s",' % (name))
        elif 'from cscripts.csworkflow    import csworkflow' in line:
            print('from cscripts.%s    import %s' % (name, name))
        print(line,end=''),

    # Update cscripts/Makefile.am file
    for line in fileinput.FileInput('cscripts/Makefile.am',inplace=1):
        if '$(srcdir)/csworkflow.py' in line:
            print('                $(srcdir)/%s.py \\' % (name))
        elif '$(srcdir)/csworkflow.par' in line:
            print('           $(srcdir)/%s.par \\' % (name))
        elif 'csworkflow \\' in line:
            print('              %s \\' % (name))
        elif '$(top_srcdir)/test/test_csworkflow.py' in line:
            print('              $(top_srcdir)/test/test_%s.py \\' % (name))
        print(line,end=''),

    # Update test/test_python_cscripts.py file
    for line in fileinput.FileInput('test/test_python_cscripts.py',inplace=1):
        if 'import test_csworkflow' in line:
            print('import test_%s' % (name))
        elif 'test_csworkflow.Test()' in line:
            print('             test_%s.Test(),' % (name))
        print(line,end=''),

    # Update reference manual index file
    for line in fileinput.FileInput('doc/source/users/reference_manual/index.rst',inplace=1):
        if 'csbkgmodel --- Generates background model for 3D analysis <csbkgmodel>' in line:
            print('   %s --- ToDo: Describe what script is doing <%s>' % (name, name))
        print(line,end=''),

    # Return
    return


# ========== #
# ctool menu #
# ========== #
def ctool_menu():
    """
    ctool menu
    """
    # Annonce actions
    print("")
    print("Add ctool")
    print("---------")

    # Stay in loop until there is a final confirmation
    while True:

        # Enter tool name
        ctoolname = response('Please enter a ctool name (e.g. "ctpntsim")')

        # Enter baseclass
        print('From which baseclass should the ctool derive?')
        print('[1] ctool')
        print('[2] ctobservation')
        print('[3] ctlikelihood')
        waiting = True
        choice  = 1
        while waiting:
            choice = str(input('Enter your choice: '))
            if choice == '1' or choice == '2' or choice == '3':
                waiting = False
        if choice == '1':
            basename  = 'base'
            baseclass = 'ctool'
        elif choice == '2':
            basename  = 'obs'
            baseclass = 'ctobservation'
        elif choice == '3':
            basename  = 'like'
            baseclass = 'ctlikelihood'

        # Enter further information
        what        = response('Please say what the tool is for (e.g. "Pointing simulation")')
        author      = response('Please enter your name (e.g. "Joe Public")')
        email       = response('Please enter your e-mail (e.g. "joe.public@dot.com")')
        affiliation = response('Please enter your affiliation (e.g. "ESA")')

        # Ask to confirm module summary
        print('\nAll right. Have now:')
        print('ctool name .......: "%s"' % ctoolname)
        print('Base class .......: "%s"' % baseclass)
        print('Tools descriptor .: "%s"' % what)
        print('Your name ........: "%s"' % author)
        print('Your e-mail ......: "%s"' % email)
        print('Your affiliation .: "%s"' % affiliation)
        if confirm('Is this correct?'):
            break

    # Set tokens
    tokens = set_ctool_tokens(ctoolname, basename, what, author, email, affiliation)

    # Add ctool
    add_ctool(ctoolname, tokens, basename)

    # Return
    return


# ============ #
# cscript menu #
# ============ #
def cscript_menu():
    """
    cscript menu
    """
    # Annonce actions
    print("")
    print("Add cscript")
    print("-----------")

    # Stay in loop until there is a final confirmation
    while True:

        # Enter tool name
        cscriptname = response('Please enter a cscript name (e.g. "cspntsim")')

        # Enter baseclass
        print('From which baseclass should the cscript derive?')
        print('[1] cscript')
        print('[2] csobservation')
        print('[3] cslikelihood')
        waiting = True
        choice  = 1
        while waiting:
            choice = str(input('Enter your choice: '))
            if choice == '1' or choice == '2' or choice == '3':
                waiting = False
        if choice == '1':
            basename  = 'base'
            baseclass = 'cscript'
        elif choice == '2':
            basename  = 'obs'
            baseclass = 'csobservation'
        elif choice == '3':
            basename  = 'like'
            baseclass = 'cslikelihood'

        # Enter further information
        what        = response('Please say what the script is for (e.g. "Pointing simulation")')
        author      = response('Please enter your name (e.g. "Joe Public")')
        email       = response('Please enter your e-mail (e.g. "joe.public@dot.com")')
        affiliation = response('Please enter your affiliation (e.g. "ESA")')

        # Ask to confirm module summary
        print('\nAll right. Have now:')
        print('cscript name .....: "%s"' % cscriptname)
        print('Base class .......: "%s"' % baseclass)
        print('Script descriptor : "%s"' % what)
        print('Your name ........: "%s"' % author)
        print('Your e-mail ......: "%s"' % email)
        print('Your affiliation .: "%s"' % affiliation)
        if confirm('Is this correct?'):
            break

    # Set tokens
    tokens = set_cscript_tokens(cscriptname, basename, what, author, email, affiliation)

    # Add cscript
    add_cscript(cscriptname, tokens, basename)

    # Return
    return


# ================ #
# Manage main menu #
# ================ #
def main_menu():
    """
    Manage main menu
    """
    # Print main menu
    print('[1] Add ctool')
    print('[2] Add cscript')
    print('[q] Quit')

    # Wait for the input
    waiting = True
    while waiting:
        choice = str(input('Enter your choice: '))
        if choice == '1' or choice == '2' or choice == 'q':
            waiting = False

    # Return choice
    return choice


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Check if script has been started from ctools root directory
    if not os.path.isfile('ctools.pc.in'):
        print('"codgen.py" script needs to be started from ctools source code '
              'root directory, you are somewhere else. Quit now.')
        sys.exit()

    # Clear console
    os.system('clear')

    # Print header
    print('ctools code generator')
    print('=====================')
    print('')

    # Enter endless loop
    while True:

        # Show main menu
        choice = main_menu()

        # Dispatch according to choice
        if choice == '1':
            ctool_menu()
            print('')
        elif choice == '2':
            cscript_menu()
            print('')
        elif choice == 'q':
            break
