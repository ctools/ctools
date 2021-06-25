.. _comobsback:

comobsback
==========

Generate background model for COMPTEL observations.


Synopsis
--------

This script generates a background model for COMPTEL observations.


General parameters
------------------

``inobs [file]``
    Input observation definition XML file.

``inmodel [file]``
    Input model definition XML file. If ``NONE`` is specified the input model
    will be ignored.

``(suffix = srclix) [string]``
    Suffix for DRB files. This suffix will be appended to each DRB file so that
    distinctive file names can be created in the same outfolder.

``(outfolder = dri) [string]``
    Output folder for DRB files.

``outobs [file]``
    Output observation definition XML file.

``bkgmethod <PHINOR|BGDLIXA|BGDLIXE> [string]``
    Method for background computation.

``(nrunav = 3) [integer]``
    Number of bins used for running average.

``(navgr = 3) [integer]``
    Number of bins used for averaging.

``(nincl = 13) [integer]``
    Number of Phibar layers to include.

``(nexcl = 0) [integer]``
    Number of Phibar layers to exclude.


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

``(logfile = comobsback.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

