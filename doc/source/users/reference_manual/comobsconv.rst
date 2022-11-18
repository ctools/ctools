.. _comobsconv:

comobsconv
==========

Convolves models with COMPTEL reponse.


Synopsis
--------

This script convolves celestial models of gamma-ray emission with the COMPTEL
response and stores the result in a FITS file. Only model component for which
all spatial parameters are fixed will be considered. The response will be
stored in FITS files using a compressed format that is also used internally by
the software to cache response computations. The convolved response information
is attached using a ``RSP`` tag to the output observation definition XML file
and will be automatically used for any further processing.

The file names of the convolved response files will be constructed from the
background model file names by replacing the string ``drb`` with the string
specified using the hidden ``filetype`` parameter. The convolved response
files will be stored in the folder specified by the hidden ``outfolder``
parameter. Once written, a convolved response file will never be written again.
If a convolved response file should be rewritten, the ``filetype`` parameter
needs to modified so that a different filename will be used.


General parameters
------------------

``inobs [file]``
    Input observation definition XML file.

``inmodel [file]``
    Input model definition XML file.

``(outfolder = dri) [string]``
    Output folder for convolved response files.

``outobs [file]``
    Output observation definition XML file.

``(filetype = rsp) [string]``
    File type.


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

``(logfile = comobsconv.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

None
