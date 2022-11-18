.. _comobssim:

comobssim
=========

Simulate COMPTEL observations.


Synopsis
--------

This script simulates DRE data for an observation based on a specified model.

If ``add=yes`` then the simulated events are added on top of an existing DRE. This
mode can be used to add for example a simulated source to some real data in
order to check whether the simulated source parameters can be recovered. In
this case, inmodel needs to specify only the source to be added.

If ``add=no`` then any preexisting events will be dropped, and the produced DRE
will only contain the simulated events. This mode is useful if also a background
model should be simulated.

The random number generator seed may be set via the ``seed`` parameter. A given
seed value will always produce a given and reproducible sequence of events. If
multiple DRE files should be produced that differ randomly in their events, each
simulation should be done using a different seed value.

The produced DRE files will be written in the folder specified by ``outfolder``,
and the tool will produce an output observation definition XML file, specified
by ``outobs``, that points to the simulated DRE files.


General parameters
------------------

``inobs [file]``
    Input observation definition XML file.

``inmodel [file]``
    Input model definition XML file.

``(suffix = ) [string]``
    Suffix for DRE files. This suffix will be appended to each DRE file so that
    distinctive file names can be created in the same outfolder.

``(outfolder = dri) [string]``
    Output folder for DRI files.

``outobs [file]``
    Output observation definition XML file.

``add [boolean]``
    Add to existing events?

``seed [integer]``
    Random number generator seed.


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

``(logfile = comobssim.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

None
