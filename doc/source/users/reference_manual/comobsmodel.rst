.. _comobsmodel:

comobsmodel
===========

Generate model for binned COMPTEL observations.


Synopsis
--------

This script generates a model for COMPTEL data analysis.


General parameters
------------------

``inobs [file]``
    Input observation definition XML file.

``outmodel [file]``
    Output model definition XML file.

``ra [real]``
    Right Ascension of point source (deg). Specify NONE if no point source
    should be added to the model.

``dec [real]``
    Declination of point source (deg). Specify NONE if no point source should
    be added to the model.

``srcname [real]``
    Name of point source.

``brems <NONE|MAP|CUBE> [string]``
    Bremsstrahlung component.

``ic <NONE|MAP|CUBE> [string]``
    Inverse Compton component.

``iso <NONE|MAP|CUBE> [string]``
    Isotropic component.

``diffusetype <NODES|BINS> [string]``
    Diffuse model type.

``bkgtype <NODES|BINS> [string]``
    Background model type.

``(bremsmap = $COMDATA/../skymaps/galprop/map_comptel_bremsstrahlung.fits) [file]``
    Bremsstrahlung map file name.

``(bremscube = $COMDATA/../skymaps/galprop/bremss_mapcube_54_77Xvarh7S.fits) [file]``
    Bremsstrahlung cube file name.

``(icmap = $COMDATA/../skymaps/galprop/map_comptel_ic.fits) [file]``
    Inverse Compton map file name.

``(iccube = $COMDATA/../skymaps/galprop/ics_isotropic_mapcube_54_77Xvarh7S.fits) [file]``
    Inverse Compton cube file name.


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

``(logfile = comobsmodel.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

