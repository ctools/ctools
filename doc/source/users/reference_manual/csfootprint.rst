.. _csfootprint:

csfootprint
===========

Generates a carbon footprint report for the usage of ctools on your computed.


Synopsis
--------

This script generates a carbon footprint report based on the usage of ctools
on your computed. ctools usage statistics are collected in the file
``~/.gamma/statistics.xml`` that is located in the user's home directory. The
information file is updated using a dedicated daemon that works in the background
and collects information about ctools usage.


General parameters
------------------

``(infile = $HOME/.gamma/statistics.xml) [file]``
    Input high-level statistics XML file.

``outfile [file]``
    Output graphics file (``NONE`` if no graphics should be generated).

``tmin [time]``
    Start time for footprint reporting (UTC string, JD, MJD or MET in seconds).
    Start times given in MET seconds are counted with respect to the time
    reference of the input observation(s).
    If ``INDEF``, ``NONE``, ``UNDEF`` or ``UNDEFINED`` is passed as value, the
    report will be generated for the full data provided in the statistics
    file.

``tmax [time]``
    Stop time for footprint reporting (UTC string, JD, MJD or MET in seconds).
    Stop times given in MET seconds are counted with respect to the time
    reference of the input observation(s).
    If ``INDEF``, ``NONE``, ``UNDEF`` or ``UNDEFINED`` is passed as value, the
    report will be generated for the full data provided in the statistics
    file.


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
    Specifies whether an existing output file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = csfootprint.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

None
