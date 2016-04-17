#! /usr/bin/env python
# ==========================================================================
# This script executes a workflow defined in an XML file.
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
import gammalib
import ctools
import cscripts
import sys


# ================ #
# csworkflow class #
# ================ #
class csworkflow(ctools.cscript):
    """
    """
    def __init__(self, *argv):
        """
        Constructor.
        """
        # Set name and version
        self.name    = "csworkflow"
        self.version = "1.0.0"

        # Set members
        self.workflow = gammalib.GXml()
        self.actors   = []

        # Make sure that parfile exists
        self.parfile()

        # Initialise application
        if len(argv) == 0:
            ctools.cscript.__init__(self, self.name, self.version)
        elif len(argv) ==1:
            ctools.cscript.__init__(self, self.name, self.version, *argv)
        else:
            raise TypeError("Invalid number of arguments given.")

        # Set logger properties
        self.log_header()
        self.log.date(True)

        # Return
        return

    def __del__(self):
        """
        Destructor.
        """
        # Return
        return

    def parfile(self):
        """
        Check if parfile exists. If parfile does not exist then create a
        default parfile. This kluge avoids shipping the cscript with a parfile.
        """
        # Set parfile name
        parfile = self.name+".par"

        # Load parfile or create one
        try:
            pars = gammalib.GApplicationPars(parfile)
        except:
            # Signal if parfile was not found
            print("Parfile "+parfile+" not found. Create default parfile.")

            # Create default parfile
            pars = gammalib.GApplicationPars()
            pars.append(gammalib.GApplicationPar("inflow","f","h","test/data/workflow.xml","","","Input workflow XML file"))
            pars.append_standard()
            pars.append(gammalib.GApplicationPar("logfile","f","h","csworkflow.log","","","Log filename"))
            pars.save(parfile)

        # Return
        return

    def get_parameters(self):
        """
        Get parameters from parfile.
        """
        # Load workflow XML file
        xmlfile       = self["inflow"].filename()
        self.workflow = gammalib.GXml(xmlfile.url())

        # Return
        return

    def execute(self):
        """
        Execute the script.
        """
        # Run the script
        self.run()

        # Return
        return

    def run(self):
        """
        Run the script.
        """
        # Switch screen logging on in debug mode
        if self.logDebug():
            self.log.cout(True)

        # Get parameters
        self.get_parameters()

        #  Write input parameters into logger
        if self.logTerse():
            self.log_parameters()
            self.log("\n")

        # Parse XML file
        if self.logTerse():
            self.log("\n")
            self.log.header1("Parse workflow XML file")
        self.parse_workflow()

        # Execute workflow
        if self.logTerse():
            self.log("\n")
            self.log.header1("Execute workflow")
        self.execute_workflow()

        # Return
        return

    def parse_workflow(self):
        """
        """
        # Get workflow
        workflow = self.workflow.element("workflow")
        
        # Get number of actors
        num_actors = workflow.elements("actor")

        # Initialise actors
        self.actors = []

        # Loop over all actors
        for i in range(num_actors):
        
            # Get actor
            actor = workflow.element("actor", i)

            # Initialise parameters and input actors
            input_parameters  = []
            output_parameters = []
            input_actors      = []
            output_actors     = []

            # Get actor attributes
            name = actor.attribute("name")
            tool = actor.attribute("tool")

            # Get actor input parameters
            if actor.elements("input") > 0:
                actor_inputs = actor.element("input")
                num_inputs   = actor_inputs.elements("parameter")
                for k in range(num_inputs):
                    input       = actor_inputs.element("parameter", k)
                    input_name  = input.attribute("name")
                    input_value = input.attribute("value")
                    input_actor = input.attribute("actor")
                    parameter   = {'name': input_name, \
                                   'value': input_value, \
                                   'actor': input_actor}
                    input_parameters.append(parameter)
                    if input_actor != "":
                        if input_actor not in input_actors:
                            input_actors.append(input_actor)

            # Get actor output parameters
            if actor.elements("output") > 0:
                actor_output = actor.element("output")
                num_outputs  = actor_output.elements("parameter")
                for k in range(num_outputs):
                    output       = actor_output.element("parameter", k)
                    output_name  = output.attribute("name")
                    output_value = output.attribute("value")
                    output_actor = output.attribute("actor")
                    parameter    = {'name': output_name, \
                                    'value': output_value, \
                                    'actor': output_actor}
                    output_parameters.append(parameter)
                    if output_actor != "":
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
            self.actors.append(entry)

            # Log information about actors
            if self.logNormal():
                self.log.parformat("Actor \""+name+"\"")
                self.log(tool)
                self.log("\n")
                if num_inputs == 0:
                    self.log.parformat("  Predecessor")
                    self.log("none\n")
                else:
                    if num_inputs == 1:
                        self.log.parformat("  Predecessor")
                    else:
                        self.log.parformat("  Predecessors")
                    for k in range(num_inputs):
                        if k > 0:
                            self.log(", ")
                        self.log("\""+input_actors[k]+"\"")
                    self.log("\n")

        # Return
        return

    def execute_workflow(self):
        """
        """
        # Continue while there are actors
        while len(self.actors) > 0:

            # Find actors which are ready
            num_executed = 0
            for actor in self.actors:

                # If actor is ready then execute it
                if actor['status'] == 'ready':

                    # Log execution start
                    if self.logNormal():
                        self.log.parformat("Execute actor")
                        self.log("\""+actor['name']+"\"\n")

                    # Execute actor
                    self.execute_actor(actor)

                    # Log execution finish
                    if self.logNormal():
                        self.log.parformat("Finished actor execution")
                        self.log("\""+actor['name']+"\"\n")

                    # Set actor status to finished
                    actor['status'] = 'finished'

                    # Increment number of executed actors
                    num_executed += 1

            # Break if no actors have been executed
            if num_executed < 1:
                break
            
            # Update actor status
            for actor in self.actors:
                if actor['status'] != 'finished':
                    input_actors = actor['input_actors']
                    ready        = True
                    for input_actor in input_actors:
                        a = self.get_actor(input_actor)
                        if a['status'] != 'finished':
                            ready = False
                            break
                    if ready:
                        actor['status'] = 'ready'

        # Return
        return

    def execute_actor(self, actor):
        """
        Excute an actor.
        """
        # Log input parameters
        if self.logNormal():
            pars = actor['input_parameters']
            for par in pars:
                self.log.parformat("  Input parameter")
                self.log(par['name'])
                self.log("=")
                self.log(self.get_parameter_value(par))
                self.log("\n")

        # Log output parameters
        if self.logNormal():
            pars = actor['output_parameters']
            for par in pars:
                self.log.parformat("  Output parameter")
                self.log(par['name'])
                self.log("=")
                self.log(self.get_parameter_value(par))
                self.log("\n")

        # Set actor tool
        if 'tool' in actor:
            tool = actor['tool']
            if tool[:2] == 'ct':
                tool_eval = 'ctools.'+tool+'()'
            elif tool[:2] == 'cs':
                tool_eval = 'cscripts.'+tool+'()'
            else:
                tool_eval = ""
        if tool_eval == "":
            return

        # Create actor tool object
        object = eval(tool_eval)

        # Set actor input parameters
        pars = actor['input_parameters']
        for par in pars:
            app_par = object[par['name']]
            app_par.value(self.get_parameter_value(par))
            object[par['name']] = app_par

        # Set actor output parameters
        pars = actor['output_parameters']
        for par in pars:
            app_par = object[par['name']]
            app_par.value(self.get_parameter_value(par))
            object[par['name']] = app_par

        # Execute actor
        object.logFileOpen()
        object.execute()
        
        # Return
        return

    def get_parameter_value(self, par):
        """
        Get the parameter value from the respective actor.
        """
        # If the actor is empty then simply return the value
        if (par['actor'] == ""):
            value = par['value']

        # ... otherwise get the value from the actor
        else:
            found = False
            actor = self.get_actor(par['actor'])
            for output in actor['output_parameters']:
                if output['name'] == par['value']:
                    value = output['value']
                    found = True
                    break
            if not found:
                msg = "Parameter \""+par['value']+"\" is not an output "\
                      "parameter of actor \""+par['actor']+"\"."
                raise RuntimeError(msg)
        
        # Return value
        return value

    def get_actor(self, name):
        """
        Get the actor with the specified name.
        """
        # Return actor with specified name
        for actor in self.actors:
            if actor['name'] == name:
                return actor

        # Throw an exception
        msg = "Actor \""+name+"\" not found."
        raise RuntimeError(msg)
 
        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Executes workflow.
    """
    # Create instance of application
    app = csworkflow(sys.argv)

    # Open logfile
    app.logFileOpen()

    # Execute application
    app.execute()
