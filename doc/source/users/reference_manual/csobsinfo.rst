.. _csobsinfo:

csobsinfo
=========

Dumps information about an observation container into log file or on 
screen.


Synopsis
--------

This script provides detailed information about a given observation container.
For IACT observations usually several observations with different configurations
are combined. The script computes e.g. the overall energy range which could be
an important input for :doc:`ctbin` to create a suitable energy range. In
addition the script prints the pointing attributes like the average
zenith/azimuth angle. If the hidden boolean parameter "offset" is specified to
"yes", the script queries for the target position and computes the run-by-run
offset as well as its average. In the observation summary, the time range of
the observations including total ontime and (dead-time corrected) livetime is
shown.

The script further contains the hidden parameter "ds9file". If specified with
a file name, the tool writes out a ds9 region file including a cross for each
pointing position in the sky.  

When executed from within python, the script provides methods to return arrays
containing observation-wise parameters. This might be useful for graphical
display within a script. The following methods are implemented:
- zeniths (returning all zenith angles)
- azimuths (returning all azimuth angles)
- offsets (returning all offset angles, provided the hidden parameter offset=yes)
- ebounds (returning energy boundaries)
- gtis (returning good time intervals)


General parameters
------------------

``inobs [file]``
    Event list, counts cube, or observation definition file

``outds9file [file]``
    DS9 region file containing pointing directions

``(offset = no) [file]``
    Compute offset from target to pointing positions

``ra [real]``
    Target right ascension for offset computation (only queried if offset=yes)

``dec [real]``
    Target declination for offset computation (only queried if offset=yes)


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
    Specifies whether an existing output file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csobsinfo.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

None
