.. _ctobssim:

ctobssim
========

Simulates CTA events.


Synopsis
--------

Generate photon events from astrophysical sources and background events 
from an instrumental background model.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file. If provided (i.e. the
    parameter is not blank or NONE), the pointing definition and eventually the
    response information will be extracted from the input file for event
    simulation.

``inmodel [file]``
    Input model XML file.
 	 	 
``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.
 	 	 
``(edisp = no) [boolean]``
    Apply energy dispersion?
 	 	 
``outevents [file]``
    Output event list or observation definition XML file.
 	 	 
``(prefix = sim_events_) [string]``
    Prefix for event list in observation definition XML file.
 	 	 
``(seed = 1) [integer]``
    Integer seed value to be used for Monte Carlo simulations. Keep this 
    parameter at the same value for repeatable simulations, or increment 
    this value for subsequent runs if non-repeatable simulations are
    required.
 	 	 
``ra [real]``
    Right Ascension of CTA pointing (J2000, in degrees).
 	 	 
``dec [real]``
    Declination of CTA pointing (J2000, in degrees).
 	 	 
``rad [real]``
    Radius of CTA field of view (simulation cone radius) (in degrees).
 	 	 
``tmin [real]``
    Mission elapsed start time of observation (in seconds).
 	 	 
``tmax [real]``
    Mission elapsed stop time of observation (in seconds).
 	 	 
``emin [real]``
    Lower energy limit of simulated events (in TeV).
 	 	 
``emax [real]``
    Upper energy limit of simulated events (in TeV).
 	 	 
``(deadc = 0.95) [real]``
    Average deadtime correction factor.

``(maxrate = 1.0e6) [real]``
    Maximum photon rate for source models. Source models that exceed this
    maximum photon rate will lead to an exception as very likely the
    specified model normalisation is too large (probably due to the
    a misinterpretation of units). Note that ctools specifies intensity
    units per MeV.


Standard parameters
-------------------

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output counts cube should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctobssim.log) [string]``
    Name of log file.


Related tools
-------------

None
