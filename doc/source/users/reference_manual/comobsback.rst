.. _comobsback:

comobsback
==========

Generate background model for COMPTEL observations.


Synopsis
--------

This script generates a background model for COMPTEL observations using either
the PHINOR, the BGDLIXA or the BGDLIXE algorithm. The resulting background
model will be written into a DRB FITS file that will be stored in the folder
defined by the ``outfolder`` parameter. The DRB FITS file name will be derived
from the DRE filename with the background modelling method and parameters
attached. An additional suffix, specified by the ``suffix`` parameter, may be
added to the DRB filename.

The script also generates an output observation definition XML file, specified
by the ``outobs`` parameter, that contains a link to the generated background
model. This output observation definition XML file should be used for further
data analysis.


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
    Number of Chi/Psi bins used for running average (relevant for ``BGDLIXA``
    method).

``(navgr = 9) [integer]``
    Number of Chi/Psi bins used for averaging (relevant for ``BGDLIXA`` and
    ``BGDLIXE`` methods).

``(nincl = 5) [integer]``
    Number of Phibar layers to include (relevant for ``BGDLIXA`` and ``BGDLIXE``
    methods).

``(nexcl = 0) [integer]``
    Number of Phibar layers to exclude (relevant for ``BGDLIXA`` and ``BGDLIXE``
    methods).

``(phinor = yes) [boolean]``
    Specifies whether the Phibar distribution of the resulting background model
    should be normalised to the number of observed events minus the number of
    source events.


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

None
