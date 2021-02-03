#! /usr/bin/env python
# ==========================================================================
# This script tests the Sphinx tutorials
#
# Usage:
#   ./tutorials.py
#
# --------------------------------------------------------------------------
#
# Copyright (C) 2021 Juergen Knoedlseder
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
import os
import sys
import subprocess
import gammalib
import ctools
import cscripts


# ===================================== #
# Test class for tutorials verification #
# ===================================== #
class tutorials(gammalib.GPythonTestSuite):
    """
    Test class for tutorials verification
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Initialise results
        self.results = None

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('Tutorials Verification')

        # Append tutorials
        self.append(self.tutorials_1dc, 'Test 1DC tutorial')

        # Return
        return

    # Test line execution
    def _execute_cmd(self, cmd, args):
        """
        Execute command
        """
        # Create subprocess for tool
        p = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE,
                                             stdin=subprocess.PIPE,
                                             stderr=subprocess.PIPE)

        # Write command arguments parameters
        if len(args) > 0:
            p.stdin.write(b'%s' % args)
            p.stdin.close()

        # Get results
        res = p.communicate()

        # Set error flag and text
        name = 'Test execution of command %s.' % cmd
        self.test_assert((len(res[1]) == 0), name, res[1])

        # Return
        return

    # Test line execution
    def _test_line_execution(self, lines):
        """
        Test the execution of a sequence of commands
 
        Parameters
        ----------
        lines : list of str
            Command lines
        """
        # Initialise commands and arguments
        cmd  = ''
        args = ''
        
        # Loop over all lines
        for line in lines:

            # If line starts with a string then extract or execute command
            if line[0:2] == '$ ':

                # If we have a command then execute it
                if cmd != '':
                    print(cmd, args)
                    self._execute_cmd(cmd, args)
                
                # Set new command
                cmd  = line[2:]
                args = ''

            # ... otherwise add argument
            else:
                args += '%s\n' % (gammalib.strip_whitespace(line))

        # If there is a pending command then execute it
        if cmd != '':
            print(cmd, args)
            self._execute_cmd(cmd, args)

        # Return
        return

    # Test ctool or cscript execution
    def _test_ctool_execution(self, tool, lines):
        """
        Test the execution of a ctool or cscript

        The methods test the execution of a ctool or cscript as executed from
        the command line emulating the entering of the queried parameters. Any
        error that occurs in the execution will be signaled and the Python
        traceback will be put in the error message field.

        Parameters
        ----------
        tool : str
            ctool name
        lines : list of str
            Command lines
        """
        # Initialise command
        cmd = ''

        # Append lines
        for line in lines:
            pos = line.find(']')
            if pos != -1:
                arg  = gammalib.strip_whitespace(line[pos+1:])
                cmd += '%s\n' % arg

        # Create subprocess for tool
        p = subprocess.Popen([tool], stdout=subprocess.PIPE,
                                     stdin=subprocess.PIPE,
                                     stderr=subprocess.PIPE)

        # Write command line parameters
        p.stdin.write(b'%s' % cmd)

        # Get results
        res = p.communicate()

        # Close pipeline
        p.stdin.close()

        # Set error flag and text
        name = 'Test execution of %s.' % tool
        self.test_assert((len(res[1]) == 0), name, res[1])

        # Return
        return

    # Add workflow step
    def _add_workflow_step(self, tool, lines):
        """
        Add workflow step

        Parameters
        ----------
        tool : str
            Tool name
        lines : list of str
            Command lines

        Returns
        -------
        cmd : str
            Workflow command
        """
        # If tool is a ctool or cscript then add a ctool or cscript workflow
        # step
        if tool[0:2] == 'ct' or tool[0:2] == 'cs':
            self._test_ctool_execution(tool, lines)

        # ... otherwise execute lines
        else:
            all_lines = ['$ %s' % tool]
            all_lines.extend(lines)
            self._test_line_execution(all_lines)

        # Return
        return

    # Extract workflow
    def _get_workflow(self, filename):
        """
        Extract workflow from Sphinx rst file

        Parameters
        ----------
        filename : str
            RST filename

        Returns
        -------
        workflow : `~gammalib.GXml`
            Workflow XML file
        """
        # Open RST file
        f = open(filename, 'r')

        # Signal code block
        code = False

        # Loop over lines
        for line in f:

            # If we are not in a code block then check whether a code block
            # starts
            if not code:
                pos = line.find('code-block:: bash')
                if pos != -1:
                    code  = True
                    start = -1
                    tool  = ''
                    lines = []
                    continue

            # ... otherwise
            else:

                # If we have no tool name then find tool name
                if tool == '':
                    pos = line.find(' $ ')
                    if pos != -1:
                        start = pos+1
                        tool  = line[pos+3:-1]
                        continue

                # ... otherwise append lines or add workflow step
                else:
                    blank = gammalib.strip_whitespace(line[0:start])
                    cmd   = line[start:-1]
                    if len(blank) == 0 and cmd[0] != ' ':
                        lines.append(cmd)
                        continue
                    else:
                        self._add_workflow_step(tool, lines)
                        code = False
                        continue

        # Add step if some code is pending
        if code:
            self._add_workflow_step(tool, lines)

        # Return
        return

    # Test 1DC tutorials
    def tutorials_1dc(self):
        """
        Test 1DC tutorials
        """
        # Set environment variables
        os.environ['CTADATA'] = '/project-data/cta/data/1dc'
        os.environ['CALDB']   = '/project-data/cta/data/1dc/caldb'

        # Extract workflow
        #workflow = self._get_workflow('../doc/source/users/tutorials/1dc/first_onoff.rst')
        workflow = self._get_workflow('first_onoff.rst')

        # Run workflow

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Allocate test suite container
    suites = gammalib.GTestSuites('ctools tutorials verification')

    # Allocate test suite and append it to the container
    suite_tutorials = tutorials()

    # Setup test suit
    suite_tutorials.set()

    # Append test suite to container
    suites.append(suite_tutorials)

    # Create pfiles directory
    try:
        os.mkdir('pfiles')
    except:
        pass

    # Copy ctools parameter files into pfiles directory
    os.system('cp -r ../src/*/*.par pfiles/')

    # Set PFILES environment variable
    os.environ['PFILES'] = 'pfiles'

    # Run test suite
    success = suites.run()

    # Save test results
    suites.save('reports/tutorials.xml')

    # Set return code
    if success:
        rc = 0
    else:
        rc = 1

    # Exit with return code
    sys.exit(rc)
