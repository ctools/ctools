#! /usr/bin/env python
# ==========================================================================
# Executes analysis workflow
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
import gammalib
import ctools
import cscripts


# ================ #
# csworkflow class #
# ================ #
class csworkflow(ctools.cscript):
    """
    Executes an analysis workflow.

    The ``csworkflow`` script executes an analysis workflow defined in an
    XML file.
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor.

        Parameters
        ----------
        argv : list of str
            List of IRAF command line parameter strings of the form
            ``parameter=3``.

        Raises
        ------
        TypeError
            An invalid number of command line arguments was provided.
        """
        # Set name and version
        self._name    = 'csworkflow'
        self._version = '1.2.0'

        # Set members
        self._workflow = gammalib.GXml()
        self._actors   = []

        # Initialise application by calling the appropriate class constructor
        self._init_cscript(argv)

        # Return
        return


    # Private methods
    #
    # Query all user parameters
    def _get_parameters(self):

        # Load workflow XML file
        xmlfile        = self['inflow'].filename()
        self._workflow = gammalib.GXml(xmlfile.url())

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return

    # Parse the workflow XML file
    def _parse_workflow(self):

        # Get workflow
        workflow = self._workflow.element('workflow')
        
        # Get number of actors
        num_actors = workflow.elements('actor')

        # Initialise actors
        self._actors = []

        # Loop over all actors
        for i in range(num_actors):
        
            # Get actor
            actor = workflow.element('actor', i)

            # Initialise parameters and input actors
            input_parameters  = []
            output_parameters = []
            input_actors      = []
            output_actors     = []

            # Get actor attributes
            name = actor.attribute('name')
            tool = actor.attribute('tool')

            # Get actor input parameters
            if actor.elements('input') > 0:
                actor_inputs = actor.element('input')
                num_inputs   = actor_inputs.elements('parameter')
                for k in range(num_inputs):
                    input       = actor_inputs.element('parameter', k)
                    input_name  = input.attribute('name')
                    input_value = input.attribute('value')
                    input_actor = input.attribute('actor')
                    parameter   = {'name': input_name, \
                                   'value': input_value, \
                                   'actor': input_actor}
                    input_parameters.append(parameter)
                    if input_actor != '':
                        if input_actor not in input_actors:
                            input_actors.append(input_actor)

            # Get actor output parameters
            if actor.elements('output') > 0:
                actor_output = actor.element('output')
                num_outputs  = actor_output.elements('parameter')
                for k in range(num_outputs):
                    output       = actor_output.element('parameter', k)
                    output_name  = output.attribute('name')
                    output_value = output.attribute('value')
                    output_actor = output.attribute('actor')
                    parameter    = {'name': output_name, \
                                    'value': output_value, \
                                    'actor': output_actor}
                    output_parameters.append(parameter)
                    if output_actor != '':
                        if output_actor not in output_actors:
                            output_actors.append(output_actor)

            # Determine number of dependencies
            num_inputs = len(input_actors)

            # Set actor status
            if num_inputs > 0:
                status = 'waiting for input'
            else:
                status = 'ready'

            # Create actor entry
            entry = {'name': name,
                     'tool': tool,
                     'input_parameters': input_parameters,
                     'input_actors': input_actors,
                     'output_parameters': output_parameters,
                     'output_actors': output_actors,
                     'status': status}

            # Append entry
            self._actors.append(entry)

            # Log information about actors
            self._log_value(gammalib.NORMAL, 'Actor "%s"' % name, tool)

            # Compute list of predecessors
            if num_inputs == 0:
                predecessors = 'none'
            else:
                predecessors = ''
                for k in range(num_inputs):
                    if k > 0:
                        predecessors += ', '
                    predecessors += '"'+input_actors[k]+'"'
 
            # Log predecessors
            self._log_value(gammalib.NORMAL,
                            gammalib.number('  Predecessor', num_inputs),
                            predecessors)

        # Return
        return

    # Execute a workflow
    def _execute_workflow(self):

        # Continue while there are actors
        while len(self._actors) > 0:

            # Find actors which are ready
            num_executed = 0
            for actor in self._actors:

                # If actor is ready then execute it
                if actor['status'] == 'ready':

                    # Log execution start
                    self._log_value(gammalib.NORMAL, 'Execute actor',
                                    '"%s"' % actor['name'])

                    # Execute actor
                    self._execute_actor(actor)

                    # Log execution finish
                    self._log_value(gammalib.NORMAL, 'Finished actor execution',
                                    '"%s"' % actor['name'])

                    # Set actor status to finished
                    actor['status'] = 'finished'

                    # Increment number of executed actors
                    num_executed += 1

            # Break if no actors have been executed
            if num_executed < 1:
                break
            
            # Update actor status
            for actor in self._actors:
                if actor['status'] != 'finished':
                    input_actors = actor['input_actors']
                    ready        = True
                    for input_actor in input_actors:
                        a = self._get_actor(input_actor)
                        if a['status'] != 'finished':
                            ready = False
                            break
                    if ready:
                        actor['status'] = 'ready'

        # Return
        return

    # Execute an actor
    def _execute_actor(self, actor):

        # Log input parameters
        pars = actor['input_parameters']
        for par in pars:
            self._log_value(gammalib.NORMAL, '  Input parameter', '%s=%s' %
                            (par['name'], self._get_parameter_value(par)))

        # Log output parameters
        pars = actor['output_parameters']
        for par in pars:
            self._log_value(gammalib.NORMAL, '  Output parameter', '%s=%s' %
                            (par['name'], self._get_parameter_value(par)))

        # Set actor tool
        if 'tool' in actor:
            tool = actor['tool']
            if tool[:2] == 'ct':
                tool_eval = 'ctools.'+tool+'()'
            elif tool[:2] == 'cs':
                tool_eval = 'cscripts.'+tool+'()'
            else:
                tool_eval = ''
        if tool_eval == '':
            return

        # Create actor tool object
        object = eval(tool_eval)

        # Set actor input parameters
        pars = actor['input_parameters']
        for par in pars:
            app_par = object[par['name']]
            app_par.value(self._get_parameter_value(par))
            object[par['name']] = app_par

        # Set actor output parameters
        pars = actor['output_parameters']
        for par in pars:
            app_par = object[par['name']]
            app_par.value(self._get_parameter_value(par))
            object[par['name']] = app_par

        # Execute actor
        object.logFileOpen()
        object.execute()
        
        # Return
        return

    # Get the parameter value from the respective actor
    def _get_parameter_value(self, par):

        # If the actor is empty then simply return the value
        if (par['actor'] == ''):
            value = par['value']

        # ... otherwise get the value from the actor
        else:
            found = False
            actor = self._get_actor(par['actor'])
            for output in actor['output_parameters']:
                if output['name'] == par['value']:
                    value = output['value']
                    found = True
                    break
            if not found:
                msg = 'Parameter "'+par['value']+'" is not an output '+ \
                      'parameter of actor "'+par['actor']+'".'
                raise RuntimeError(msg)
        
        # Return value
        return value

    # Get the actor with the specified name
    def _get_actor(self, name):

        # Return actor with specified name
        for actor in self._actors:
            if actor['name'] == name:
                return actor

        # Throw an exception
        msg = 'Actor "'+name+'" not found.'
        raise RuntimeError(msg)
 
        # Return
        return


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

        # Parse XML file
        if self._logTerse():
            self._log('\n')
            self._log.header1('Parse workflow XML file')
        self._parse_workflow()

        # Execute workflow
        if self._logTerse():
            self._log('\n')
            self._log.header1('Execute workflow')
        self._execute_workflow()

        # Return
        return

    def execute(self):
        """
        Execute the script
        """
        # Open logfile
        self.logFileOpen()

        # Run the script
        self.run()

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csworkflow(sys.argv)

    # Execute application
    app.execute()
