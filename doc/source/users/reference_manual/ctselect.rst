.. _ctselect:

ctselect
========

Select events from event list(s).


Synopsis
--------

This tool selects events from one or several event lists. Event selection 
is based on a circular acceptance region, a time interval and an energy 
interval. In addition, any expression following the cfitsio syntax can be 
used for event selection.

Optionally, ctselect may also apply energy thresholds. If ``usethres=DEFAULT``
is specified, ctselect will extract any save thresholds from the instrument
response functions, and if they exist, will apply them to the respective 
event list. Alternatively, if ``usethres=USER`` is specified, ctselect will
extract any user thresholds from the observation definition file (attributes
``emin`` and ``emax``), and if either of the attributes exists, will apply
them to the respective event list.

If an event list is provided on input, ctselect creates a new FITS file on 
output that contains only the selected events. In case that an observation 
definition file is specified on input, ctselect creates for each event file
referenced in the observation definition file a new FITS file with the value
of ``prefix`` prepended to the file name that contains only the selected
events. In addition, a new observation definition file will be written 
that references the new FITS files.

For each event file, the event selection parameters will be writted as data
selection keywords to the FITS header. These keywords are mandatory for any
unbinned maximum likelihood analysis of the event data.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file

``outobs [file]``
    Output event list or observation definition XML file.

``(prefix = "selected_") [string]``
    Prefix for output event lists in observation definition file.

``(usepnt = no) [boolean]``
    Use pointing instead of RA/DEC parameters?

``ra [real]``
    Right Ascension of acceptance cone (or ROI) centre (J2000, in degrees).
    If ``INDEF``, ``NONE``, ``UNDEF`` or ``UNDEFINED`` is passed as value, no ROI
    selection will be performed.

``dec [real]``
    Declination of acceptance cone (or ROI) centre (J2000, in degrees).
    If ``INDEF``, ``NONE``, ``UNDEF`` or ``UNDEFINED`` is passed as value, no ROI
    selection will be performed.

``rad [real]``
    Radius of acceptance cone (or ROI) centre (in degrees).
    If ``INDEF``, ``NONE``, ``UNDEF`` or ``UNDEFINED`` is passed as value, no ROI
    selection will be performed.

``tmin [time]``
    Start time for event selection (UTC string, JD, MJD or MET in seconds).
    Start times given in MET seconds are counted with respect to the time
    reference of the input observation(s).
    If ``INDEF``, ``NONE``, ``UNDEF`` or ``UNDEFINED`` is passed as value, no time
    selection will be performed.

``tmax [time]``
    Stop time for event selection (UTC string, JD, MJD or MET in seconds).
    Stop times given in MET seconds are counted with respect to the time
    reference of the input observation(s).
    If ``INDEF``, ``NONE``, ``UNDEF`` or ``UNDEFINED`` is passed as value, no time
    selection will be performed.

``emin [real]``
    Lower energy limit of events (in TeV).
    If ``INDEF``, ``NONE``, ``UNDEF`` or ``UNDEFINED`` is passed as value, no energy
    selection will be performed.

``emax [real]``
    Upper energy limit of events (in TeV).
    If ``INDEF``, ``NONE``, ``UNDEF`` or ``UNDEFINED`` is passed as value, no energy
    selection will be performed.

``(phase = NONE) [string]``
    String to apply a phase selection. The string must contain the boundaries 
    of the phase interval to be selected separated be a colon. More than one
    interval can be specified at the same time. In this case intervals must be 
    separated by a comma. Examples of valid strings are: ``phase = 0.3:0.6``,
    ``phase = 0.3:0.6,0.8:0.9`` or ``phase = 0.8:0.2``. In the last case, events
    with phases in the intervals [0.8,1.0] and [0.0,0.2] are selected. If NONE
    is passed as value, no phase selection will be performed.

``(expr = "") [string]``
    Additional event selection expression (cfitsio syntax).

``(usethres = NONE) [string]``
    Energy threshold type (one of NONE, DEFAULT or USER).


Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the event list(s) should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether existing output files should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = ctlike.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

None
