	.. _ctprob:

ctprob
========

Compute event probability for a given model.


Synopsis
--------

This tool takes an event list file (or an XML observation file) and a source 
model as input. For each event, it evalutes the probability for that event 
to come from each of the source in the model. This probability is obtained 
evaluating the differential expected counts for each source at the event 
direction and energy and normalizing it to the total differential expected 
counts for the given model. If only one source is provided in the model, 
the probability is set to 1 for each photon.
A column for each source in the model is added to the event list in order to 
store these values. The name of each column is made by the source name with 
the prefix "PROB_".


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file
 	 	 
``outobs [file]``
    Output event list or observation definition XML file.
 	 	 
``(prefix = "prob_") [string]``
    Prefix for output event lists in observation definition file.
 	 	 
``inmodel [string]``
    Input model XML file.

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrument response function.
 	 	 
``(edisp = no) [boolean]``
    Apply energy dispersion to response computation.


Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the event list(s) should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether existing output files should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctlike.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctphase`
:doc:`csphasecrv`

