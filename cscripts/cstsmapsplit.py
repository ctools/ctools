#!/usr/bin/env python
# ==========================================================================
# Create commands to split TS map computation
#
# Copyright (C) 2016-2019 Michael Mayer
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
import math
import gammalib
import ctools


# ================== #
# cstsmapsplit class #
# ================== #
class cstsmapsplit(ctools.cslikelihood):
    """
    Generates commands to split TS map computation

    This class implements the creation of file containing commands to split
    the computation of a TS map into several machines. It derives from
    the ctools.cscript class which provides support for parameter files,
    command line arguments, and logging. In that way the Python script
    behaves just as a regular ctool. 
    """

    # Consructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Initialise application by calling the appropriate class constructor
        self._init_cslikelihood(self.__class__.__name__, ctools.__version__, argv)

        # Set data members
        self._outmap       = gammalib.GFilename()
        self._bins_per_job = 0
        self._compute_null = False
        self._outfile      = gammalib.GFilename()
        self._map          = gammalib.GSkyMap()
        self._cmd          = []
        self._srcname      = ''

        # Return
        return


    # Private methods
    def _get_parameters(self):
        """
        Get parameters from parfile and setup the observation
        """
        # If there are no observations in container then get some ...
        if self.obs().size() == 0:
            self.obs(self._get_observations())

        # Set models if we have none
        if self.obs().models().size() == 0:
            self.obs().models(self['inmodel'].filename())

        # Get TS map paramters
        self._map     = self._create_map(self.obs())
        self._srcname = self['srcname'].string()
        self._outmap  = self['outmap'].filename()

        # Get additional parameters
        self._bins_per_job = self['bins_per_job'].integer()
        self._compute_null = self['compute_null'].boolean()
        self._outfile      = self['outfile'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    def _compute_null_hypothesis(self):
        """
        Compute null hypothesis

        Returns
        -------
        logl : float
            Log-likelihood of null hypothesis
        """
        # Store original models  
        models_orig = self.obs().models()

        # Get models instance
        models0 = self.obs().models()

        # Remove test source   
        models0.remove(self._srcname)

        # Set models for null hypothesis 
        self.obs().models(models0)

        # Optimise observation container
        self.obs().optimize(self.opt())

        # Set original models to container
        self.obs().models(models_orig)

        # Return optimised log-likelihood of null hypothesis
        return -(self.opt().value())


    # Public methods
    def run(self):
        """
        Run the script
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()

        # Write information into logger
        self._log_header1(gammalib.TERSE, 'Test source')
        self._log_string(gammalib.TERSE, str(self.obs().models()[self._srcname]))

        # Set log-likelihood to zero
        logL0 = 0.0

        # Pre-compute null hypothesis if requested
        if self._compute_null:

            # Write information into logger
            self._log_header1(gammalib.TERSE, 'Compute null hypothesis')

            # Compute null hypothesis
            logL0 = self._compute_null_hypothesis()

            # Write likelihood into logger
            self._log_value(gammalib.TERSE, 'Source removed', self._srcname)
            self._log_value(gammalib.TERSE, 'Log-likelihood', repr(logL0))

        # Get parameters of TS map ctool
        pars = gammalib.GApplicationPars('cttsmap.par')

        # Compute total number of jobs
        nbins = self._map.npix()
        njobs = int(math.ceil(float(nbins) / float(self._bins_per_job)) + 0.1)   

        # Set parameters to be skipped now, we will deal with them later
        skip_pars = ['binmin', 'binmax', 'logL0', 'outmap', 'logfile']

        # Set tool name for computation
        base_command = 'cttsmap'

        # Write information into logger
        self._log_header1(gammalib.TERSE, 'Create commands')
        self._log_value(gammalib.TERSE, 'Number of cttsmap calls', njobs)

        # Loop over TS map parameters
        for par in pars:

            # Skip if we deal with them later
            if par.name() in skip_pars:
                continue

            # Skip if parameter was not queried
            if not self[par.name()].was_queried():
                continue

            # Set TS map parameter according to input from
            # this script
            par.value(self[par.name()].value())

            # Append command to set parameter
            base_command += ' ' + par.name() + '=' + par.value()

        # Append null hypothesis parameter
        base_command += ' logL0=' + repr(logL0)  

        # Set binning to start from zero
        binmin = 0
        binmax = 0

        # Clear command sequence
        self._cmd = []

        # Loop over jobs and create commands
        for job in range(njobs):

            # Set specific outmap file name
            outmap = self._outmap.url().replace('.fits', '_'+str(job)+'.fits')

            # Set bin numbers to be computed in this job
            binmin = binmax
            binmax = binmin + self._bins_per_job
            if binmax > nbins:
                binmax = nbins

            # Setup sliced command
            sliced_command = '%s binmin=%s binmax=%s outmap=%s logfile=%s' % \
                             (base_command, str(binmin), str(binmax),
                              outmap, outmap.replace('.fits','.log'))

            # If running in background is requested then append a &
            if self['run_in_bkg'].boolean():
                sliced_command += ' &'

            # Append command to list of commands
            self._cmd.append(sliced_command)

        # Write information into logger
        if self._logExplicit():
            self._log('\n')
            self._log.header2('Commands')
            for cmd in self._cmd:
                self._log(cmd)
                self._log('\n')
            self._log('\n')

        # Return
        return

    def save(self):
        """
        Save commands to ASCII file
        """
        # Get filename
        filename = self._outfile.url()

        # Open file
        f = open(filename, 'w')

        # Write commands to file
        for cmd in self._cmd:
            f.write(cmd + '\n')

        # Close file
        f.close()

        # Make file executable
        os.system('chmod +x %s' % filename)

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = cstsmapsplit(sys.argv)

    # Execute application
    app.execute()
