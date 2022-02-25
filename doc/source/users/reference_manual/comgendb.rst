.. _comgendb:

comgendb
========

Generates COMPTEL database.

Synopsis
--------

This script generates a COMPTEL database from the COMPTEL HEASARC archive. The
COMPTEL database will then subsequently be used by ``comobsselect``.

On input the script expects a directory as downloaded from the HEASARC archive
containing subdirectories named ``phase01``, ``phase02`` etc. containing the
individual viewing periods for each CGRO mission phase.

On output the script will create a COMPTEL database directory that will contains
the database FITS file ``dbase.fits`` as well as subdirectories names ``tim``
and ``xml`` containing Good Time Interval FITS files and observation definition
XML files for each COMPTEL viewing period, respectively.

Optionally a reference database in form of an ASCII file can be specified to
extract status information for each viewing period.

Furthermore, if ``download=yes`` is specified the script will automatically
download the COMPTEL HEASARC archive and store the data into the ``archive``
directory. Data files that do already exists in the ``archive`` directory are
not downloaded again.


General parameters
------------------

``archive [string]``
    COMPTEL HEASARC archive directory. This directory contains the archive as
    downloaded from https://heasarc.gsfc.nasa.gov/FTP/compton/data/comptel/.

``dbase [string]``
    COMPTEL database directory. This directory will contain the COMPTEL database
    on output.

``(refdata = NONE) [file]``
    COMPTEL reference database to extract status information.

``(quality = 200) [integer]``
    Target quality flag of data.

``(download = no) [boolean]``
    Download COMPTEL HEASARC archive.


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

``(logfile = comgendb.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:ref:`comobsselect`
