.. _csviscube:

csviscube
=========

Computes visibility cube for a given IACT array and time period.


Synopsis
--------

This script computes a visibility cube which is a sky cube spanned by Right
Ascension, Declination and zenith angle that provides the number of hours
that a given position on the sky is observable under a given zenith angle
by an IACT array during a given time period.

For the moment, the script takes into account:

- a Sun constraint, specified as a minimum zenith angle of the Sun (hidden
  parameter ``sunzenith``). If the Sun is higher in the sky than this zenith
  angle the corresponding times are excluded
- a Moon constraint, specified as a minimum zenith angle of the Moon (hidden
  parameter ``moonzenith``). If the Moon is higher in the sky than this zenith
  angle the corresponding times are excluded
- a constraint on the lunar illumination (hidden parameter ``maxfli``). If the
  Moon has an illumination that is equal or larger than this fraction the
  corresponding times are excluded

On output, the script will provide a FITS file containing the visibility cube.
In addition, the FITS file will also contain a ``VISIBILITY`` extension in
form of a binary table that provides a number of parameters, such as the
positions of the Sun and the Moon, the fraction of lunar illumation and the
Dark Time as function of time.

Use the example script ``show_viscube.py`` to display a visibility cube.


General parameters
------------------

``tmin [time]``
    Start time of observation (UTC string, JD, MJD or MET in seconds).

``tmax [time]``
    Stop time of observation (UTC string, JD, MJD or MET in seconds).

``(mjdref = 51544.5) [real]``
    Reference Modified Julian Day (MJD) for MET.

``geolon [real]``
    Geographic longitude of the IACT array in degrees.

``geolat [real]``
    Geographic latitude of the IACT array in degrees.

``(sunzenith = 105.0) [real]``
    Minimum Sun zenith angle for observations in degrees.

    If the Sun is higher in the sky than the specified zenith angle, the
    corresponding times are excluded from the visibility cube.

``(moonzenith = 90.0) [real]``
    Minimum Moon zenith angle for observations in degrees.

    If the Moon is higher in the sky than the specified zenith angle, the
    corresponding times are excluded from the visibility cube.

``(maxfli = 0.4) [real]``
    Maximum fraction of lunar illumination.

    If the Moon has an illumination that is equal or larger than this fraction,
    the corresponding times are excluded from the visibility cube.

``outfile [file]``
    Output visibility cube file.

``(binsz = 1.0) [real]``
    Visibility cube pixel size in degrees/pixel.

``(dz = 1.0) [real]``
    Visibility cube zenith angle bin size in degrees.

``(zmax = 60.0) [real]``
    Maximum zenith angle of the visibility cube in degrees.


Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the visibility cube should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing visibility cube file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = csviscube.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

None
