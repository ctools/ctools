.. _comobsselect:

comobsselect
============

Select observations from COMPTEL database.


Synopsis
--------

This script selects COMPTEL observations from the database.


General parameters
------------------

``(dbase = $COMDATA/dbase/dbase.fits) [file]``
    COMPTEL database FITS file.

``outobs [file]``
    Output observation definition XML file.

``pntselect <CIRCLE|BOX> [string]``
    Pointing selection region shape. This may either be a ``CIRCLE`` or a ``BOX``.
    Only observations with pointing directions falling into the circle or box
    are selected.

``coordsys <CEL|GAL> [string]``
    Coordinate system for the pointing selection region centre (``CEL`` - celestial,
    ``GAL`` - galactic).

``ra [real]``
    If ``coordsys=CEL``: Right Ascension of pointing selection region centre (J2000, in deg).

``dec [real]``
    If ``coordsys=CEL``: Declination of pointing selection region centre (J2000, in deg).

``glon [real]``
    If ``coordsys=GAL``: Galactic longitude of pointing selection region centre (deg).

``glat [real]``
    If ``coordsys=GAL``: Galactic latitude of pointing selection region centre (deg).

``rad [real]``
    If ``pntselect=CIRCLE``: Radius of pointing selection region circle (deg).

``width [real]``
    If ``pntselect=BOX``: Width of pointing selection region box (deg).

``height [real]``
    If ``pntselect=BOX``: Height of pointing selection region box (deg).

``tmin [time]``
    Start time for observation selection. Set to ``NONE`` if no temporal
    selection should be applied.

``tmax [time]``
    Stop time for observation selection. Set to ``NONE`` if no temporal
    selection should be applied.


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

``(logfile = comobsselect.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:ref:`comgendb`

