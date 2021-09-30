.. _csadd2caldb:

csadd2caldb
===========

Adds CTA response function to calibration database.


Synopsis
--------

This script adds a CTA response function to the calibration database. The
CTA response functions can be downloaded from CTAO Zenodo directory.


General parameters
------------------

``indir [file]``
    Input IRF folder.

``outdir [file]``
    Output caldb folder.


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

``(logfile = csadd2caldb.log) [filename]``
    Log filename.


Related tools or scripts
------------------------
