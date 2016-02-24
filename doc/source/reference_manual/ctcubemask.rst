.. _ctcubemask:

ctcubemask
==========

Mask out bins of a counts cube.


Synopsis
--------

This tool masks out specific regions from a counts cube by setting the
corresponding bin values to -1. Bins with negative values will be ignored
in a maximum likelihood analysis. A mask is applied spatially and spectrally.

ctcubemask applies a spatial mask that is comprised of a circular selection
region and a list of circular exclusion regions. The circular selection region
is specified by the sky direction of the centre and the radius of the circle.
The circular exclusion regions are specified by an ASCII file in ds9 format.
The ASCII file contains one row per exclusion region, given in the format

``circle(83.63, 21.5, 0.4)``

where ``83.63`` and ``21.5`` are the Right Ascension and Declination of 
the centre and ``0.4`` is the radius (in degrees) of the exclusion circle.

ctcubemask also selects only counts cube energy layers that are fully comprised
in the energy interval specified by the ``emin`` and ``emax`` parameters.

ctcubemask generates a counts cube FITS file that is a copy of the input 
counts cube and for which all masked bins have been set to -1.


General parameters
------------------

``inobs [file]``
    Input counts cube or observation definition XML file.

``outcube [file]``
    Output counts cube or observation definition XML file.

``(prefix = filtered_) [string]``
    Prefix for output counts cube files in case that an observation
    definition XML file has been specified on input.

``(usepnt = no) [boolean]``
    Use CTA pointing direction instead of ra/dec parameters?
 	 	 
``ra [real]``
    Right Ascension of circular selection region centre (J2000, in degrees).
 	 	 
``dec [real]``
    Declination of circular selection region centre (J2000, in degrees).

``rad [real]``
    Radius of circular selection region (in degrees).

``regfile [file]``
    Input exclusion region file in ds9 format.

``emin [real]``
    Lower energy limit (in TeV).
 	 	 
``emax [real]``
    Upper energy limit (in TeV).
 	 	 

Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the counts cube should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output counts cube file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctcubemask.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctbin`
