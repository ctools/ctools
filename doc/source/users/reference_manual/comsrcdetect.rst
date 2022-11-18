.. _comsrcdetect:

comsrcdetect
============

Detect source in TS map.


Synopsis
--------

This script detects sources in a TS map.


General parameters
------------------

``inmap [file]``
    Input TS map file.

``outmodel [file]``
    Output model definition file.

``outmap [file]``
    Output TS map for diagnostics

``(outds9file = srcdetect.reg) [file]``
    Output DS9 region file.

``threshold [real]``
    Detection TS threshold.

``(maxsrcs = 20 [integer])``
    Maximum number of sources.

``(prefix = Source [string])``
    Source name prefix.


Standard parameters
-------------------

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing energy boundaries output file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = comsrcdetect.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

None
