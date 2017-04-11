.. _ctprob:

ctprob
========

Compute event probability for a given model.


Synopsis
--------

This tool takes an event list file (or an XML observation definition file) and
a model as input and computes for each event the probability that it originates
from a specific model component. The probability is computed by by evaluating
the differential event probability for a given model and normalizing these
probabilities so that the sum of the probabilities for all model components is
unity for each event.

The tool appends for each model component a single precision column to the event
list that contains the event probabilities. The name of the column is build from
the model name prefixed with ``PROB_``.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file
 	 	 
``outobs [file]``
    Output event list or observation definition XML file.
 	 	 
``(prefix = "prob_") [string]``
    Prefix for output event lists in observation definition file.
 	 	 
``inmodel [string]``
    Input model definition XML file.

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrument response function.
 	 	 
``(edisp = no) [boolean]``
    Applies energy dispersion to response computation.


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

``(logfile = ctprob.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctphase`
:doc:`csphasecrv`

