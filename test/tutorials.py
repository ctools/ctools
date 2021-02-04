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

        # Initialise private members
        self._rst_step = 0     # Step counter
        self._rst_path = ''    # Path where the tested RST file resides

        # Return
        return

    # Clean pfiles directory
    def _clean_pfiles(self):
        """
        Clean pfiles directory
        """
        # Remove all files from pfiles directory
        os.system('rm -rf %s/*.par' % os.environ['PFILES'])

        # Return
        return

    # Set result directory
    def _set_result_dir(self, path):
        """
        Set result directory

        Parameters
        ----------
        path : str
            Result directory path
        """
        # Create result directory
        try:
            os.makedirs(path)
        except:
            pass

        # Clear result directory
        os.system('rm -rf %s/*' % path)

        # Step in result directory
        os.chdir(path)

        # Return
        return

    # Test command execution
    def _test_execute_cmd(self, cmd, args):
        """
        Test command execution

        Tests the execution of a command. Currently supported commands are:
        - cd
        - mkdir
        - nano
        - ./xxxxxx (any command xxxxxx that resides in the directory of the RST file)

        Parameters
        ----------
        cmd : str
            Command line
        args : list of str
            Additional lines following the command line
        """
        # Split command line
        cmdline = cmd.split(' ')

        # Check command line
        self.test_assert((len(cmdline) > 0), 'Test command.', cmd)

        # Continue only if valid
        if len(cmdline) > 0:

            # Extract command
            command = gammalib.strip_whitespace(cmdline[0])

            # Handle commands
            if command == 'mkdir':
                self.test_value(len(cmdline), 2, 'Test command "mkdir".')
                self.test_value(len(args), 0, 'Test that no lines follow "mkdir".')
                if len(cmdline) > 1:
                    dirname = gammalib.strip_whitespace(cmdline[1])
                    os.system('mkdir -p %s' % dirname)

            elif command == 'cd':
                self.test_value(len(cmdline), 2, 'Test command "cd".')
                self.test_value(len(args), 0, 'Test that no lines follow "cd".')
                if len(cmdline) > 1:
                    dirname = gammalib.strip_whitespace(cmdline[1])
                    os.chdir(dirname)

            elif command == 'nano':
                self.test_value(len(cmdline), 2, 'Test command "nano".')
                self.test_assert((len(args) > 0), 'Test that lines follow "nano".')
                if len(cmdline) > 1:
                    fname = gammalib.strip_whitespace(cmdline[1])
                    f = open(fname, 'wb')
                    for arg in args:
                        f.write(arg)
                    f.close()

            elif command[0:2] == './':
                scriptname = command[2:]
                self.test_assert((len(args) == 0), 'Test script "%s".' % scriptname)
                tool = '%s/%s' % (self._rst_path, scriptname)

                # Debug: show tool
                #print tool.split(' '),

                # Create subprocess for tool
                p = subprocess.Popen(tool.split(' '), stdout=subprocess.PIPE,
                                                      stdin=subprocess.PIPE,
                                                      stderr=subprocess.PIPE)

                # Get results
                res = p.communicate()

                # Set error flag and text
                name = 'Test execution of %s.' % tool
                self.test_assert((len(res[1]) == 0), name, res[1])

            else:
                name = 'Test command "%s".'   % command
                msg  = 'Unknown command "%s"' % command
                self.test_assert(False, name, msg)

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
                    self._test_execute_cmd(cmd, args)
                
                # Set new command
                cmd  = line[2:]
                args = ''

            # ... otherwise add argument
            else:
                args += '%s\n' % (gammalib.strip_whitespace(line))

        # If there is a pending command then execute it
        if cmd != '':
            self._test_execute_cmd(cmd, args)

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
        # Increment step counter
        self._rst_step += 1
        
        # Initialise arguments
        args = ''

        # Set command
        cmd = tool.split(' ')

        # Append lines
        for line in lines:
            pos = line.find(']')
            if pos != -1:
                arg  = gammalib.strip_whitespace(line[pos+1:])
                args += '%s\n' % arg

        # Debug: show tool
        #print cmd,

        # Create subprocess for tool
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   stdin=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        # Write command line parameters
        p.stdin.write(b'%s' % args)

        # Get results
        res = p.communicate()

        # Close pipeline
        p.stdin.close()

        # Set error flag and text
        name = 'Test execution of %s.' % tool
        self.test_assert((len(res[1]) == 0), name, res[1])

        # Preprend step counter to logfile
        logfile = gammalib.GFilename('%s.log' % (cmd[0]))
        if logfile.exists():
            src = logfile.url()
            dst = '%2.2d_%s' % (self._rst_step, src)
            os.system('mv %s %s' % (src, dst))

        # Return
        return

    # Test workflow step
    def _test_step(self, tool, lines):
        """
        Test workflow step

        Parameters
        ----------
        tool : str
            Tool name
        lines : list of str
            Command lines
        """
        # If tool is a ctool or cscript then test ctool or cscript execution
        if tool[0:2] == 'ct' or tool[0:2] == 'cs':
            self._test_ctool_execution(tool, lines)

        # ... otherwise test line execution
        else:
            all_lines = ['$ %s' % tool]
            all_lines.extend(lines)
            self._test_line_execution(all_lines)

        # Return
        return

    # Copy file
    def _copy_file(self, filename):
        """
        Copy file from rst directory into working directory
        """
        # Set source filename
        src = self._rst_path + '/' + filename

        # Copy file
        os.system('cp %s %s' % (src, filename))

        # Return
        return

    # Test Sphinx rst file
    def _test_rst_file(self, filename):
        """
        Test Sphinx rst file

        Parameters
        ----------
        filename : str
            Sphinx rst filename
        """
        # Set RST path
        self._rst_path = os.path.abspath(os.path.dirname(filename))

        # Open RST file
        f = open(filename, 'r')

        # Signal code block
        code = False

        # Loop over lines
        for line in f:

            # If we are not in a code block then check whether a code block
            # starts
            if not code:

                # Check for start of a code block
                pos = line.find('code-block:: bash')
                if pos != -1:
                    code  = True
                    start = -1
                    tool  = ''
                    lines = []
                    continue

                # ... otherwise check for literalinclude, and if found, copy
                # relevant file
                pos = line.find('literalinclude::')
                if pos != -1:
                    filename = gammalib.strip_whitespace(line[pos+16:-1])
                    self._copy_file(filename)
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

                # ... otherwise append lines or test step
                else:
                    blank = gammalib.strip_whitespace(line[0:start])
                    cmd   = line[start:-1]
                    if len(blank) == 0 and cmd[0] != ' ':
                        lines.append(cmd)
                        continue
                    else:
                        self._test_step(tool, lines)
                        code = False
                        continue

        # Test step if some code is pending
        if code:
            self._test_step(tool, lines)

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

    # Test 1DC tutorials
    def tutorials_1dc(self):
        """
        Test 1DC tutorials
        """
        # Continue only if CTADATA1DC environment variable is set
        if 'CTADATA1DC' in os.environ:
        
            # Set environment variables
            os.environ['CTADATA'] = os.environ['CTADATA1DC']
            os.environ['CALDB']   = os.environ['CTADATA1DC']+'/caldb'

            # Set Sphinx rst file path
            path = os.path.abspath('../doc/source/users/tutorials/1dc')

            # Set result directory
            self._set_result_dir('tutorials/1dc')

            # Get current working directory
            cwd = os.getcwd()

            # Clean pfiles
            self._clean_pfiles()

            # Initialise step counter
            self._rst_step = 0

            # Test Sphinx rst files
            self._test_rst_file('%s/first_select_obs.rst' % path)
            self._test_rst_file('%s/first_select_events.rst' % path)
            self._test_rst_file('%s/first_skymap.rst' % path)
            self._test_rst_file('%s/first_detect.rst' % path)
            self._test_rst_file('%s/first_stacked.rst' % path)
            self._test_rst_file('%s/first_fitting.rst' % path)
            self._test_rst_file('%s/first_improving.rst' % path)
            self._test_rst_file('%s/first_unbinned.rst' % path)

            # Reset working directory
            os.chdir(cwd)

            # Clean pfiles
            self._clean_pfiles()

            # Initialise step counter
            self._rst_step = 0

            # Test Sphinx rst file for On/Off analysis
            self._test_rst_file('%s/first_onoff.rst' % path)

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

    # Create tutorials/pfiles directory
    try:
        os.makedirs('tutorials/pfiles')
    except:
        pass

    # Set PFILES environment variable
    os.environ['PFILES'] = os.path.abspath('tutorials/pfiles')

    # Get current working directory
    cwd = os.getcwd()

    # Run test suite
    success = suites.run()

    # Reset working directory
    os.chdir(cwd)

    # Save test results
    suites.save('reports/tutorials.xml')

    # Set return code
    if success:
        rc = 0
    else:
        rc = 1

    # Exit with return code
    sys.exit(rc)
