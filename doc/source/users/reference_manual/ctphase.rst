.. _ctphase:

ctphase
=======

Computes the phase of each event using a temporal phase curve model.


Synopsis
--------

This tool computes for each event the phase value based on a temporal phase
curve model. The phase computation is required to derive for example the phase
curve of a pulsar or a gamma-ray binary.

The tool takes on input an event list or an observation definiton XML file and
appends on output to each event file a ``PHASE`` column containing the phase
value.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.
 	 	 
``outobs [file]``
    Output event list or observation definition XML file.
 	 	 
``(prefix = "phased_") [string]``
    Prefix for output event lists in observation definition XML file.
 	 	 
``inmodel [string]``
    Input model definition XML file. If ``NONE`` is specified the phase
    computation will be based on the ``mjd``, ``phase``, ``f0``, ``f1``, and ``f2``
    parameters.

``srcname [string]``
    Name of the source in the model definition XML file which should be used
    to compute the event phases.

``mjd [real]``
    Reference time in Modified Julian Days for phase computation (in days).

``phase [real]``
    Phase value at reference time.

``f0 [real]``
    Frequency at reference time (in Hz).

``f1 [real]``
    First frequency derivative at reference time (in s^-2).

``f2 [real]``
    Second frequency derivative at reference time (in s^-3).


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

``(logfile = ctphase.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctprob`
:doc:`csphasecrv`
