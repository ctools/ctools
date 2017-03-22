.. _csebins:

csebins
=======

Generates energy boundaries for stacked analysis.


Synopsis
--------

This script generates the energy bounaries for a stacked analysis by inspecting
the variation in the effective area and the background template of the response
functions. The script takes either an observation definition XML file or a
specific instrument response function on input.

csebins writes an energy boundary extension into a FITS file.


General parameters
------------------

``inobs [file]``
    Input observation definition XML file.

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.

``outfile [file]``
    Name of the energy boundary output file.

``emin [real]``
    Lower energy limit of energy boundaries (in TeV).
 	 	 
``emax [real]``
    Upper energy limit of energy boundaries (in TeV).

``(aeffthres = 0.2) [real]``
    Fractional change in effective area that leads to insertion of a new energy
    boundary.

``(bkgthres = 0.5) [real]``
    Fractional change in background rate that leads to insertion of a new energy
    boundary.


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
    Specifies whether an existing energy boundaries output file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csebins.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctbin`
:doc:`ctexpcube`
:doc:`ctpsfcube`
:doc:`ctedispcube`
:doc:`ctbkgcube`
