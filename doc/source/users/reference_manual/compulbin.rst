.. _compulbin:

compulbin
=========

Generates pulse profiles for a specified gamma-ray pulsar.

Synopsis
--------

This script generates pulse profiles for a specified gamma-ray pulsar by binning
the events from an EVP file. On input the script takes an observation definition
XML file and produces on output a FITS file that contains the number of events
per phase bin. Events are selected for an ARM window, specified by the ``armmin``
and ``armmax`` parameters. The position of the pulsar is extracted from the
ephemerides file.


General parameters
------------------

``inobs [file]``
    Input observation definition XML file.

``outfile [file]``
    Output pulsar phase bins file.

``psrname [string]``
    Pulsar name, as specified in the pulsar ephemerides file.

``ephemerides [file]``
    Pulsar ephemerides file. Various formats are supported: psrtime, INTEGRAL,
    Fermi and tempo2 par files.

``emin [real]``
    Minimum energy (MeV).

``emax [real]``
    Maximum energy (MeV).

``pnumbins [integer]``
    Number of phase bins.

``armmin [real]``
    Minimum ARM (deg).

``armmax [real]``
    Maximum ARM (deg).

``phimin [real]``
    Minimum Phibar angle (deg).

``phimax [real]``
    Maximum Phibar angle (deg).

``psdmin [real]``
    Minimum PSD value.

``psdmax [real]``
    Maximum PSD value.

``(zetamin = 5.0) [real]``
    Minimum Earth horizon - Phibar (zeta) angle (deg).

``(fpmtflag = 0) [integer]``
    Handling of D2 modules with failed PMT flag (0: exclude, 1: include, 2: exclude PMT).

``(d1use = 1111111) [string]``
    D1 module usage (1: use, 0: don't use).

``(d2use = 11111111111111) [string]``
    D2 module usage (1: use, 0: don't use).


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

``(logfile = compulbin.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

