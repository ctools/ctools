.. _csfindobs:

csfindobs
=========

Creates a list of observation IDs that fulfils certain criteria, e.g.
pointing close to specific source.


Synopsis
--------

This script can be used to create lists of observation IDs that share certain
properties. For instance, all observations pointing around a specific sky
position can be selected. In addition, the script is flexible to use expressions
to further constrain criteria to select observations. The user can e.g. select
all observations within a certain zenith angle range by specifying the hidden
parameter expression. For example

``csfindobs expression="ZEN_PNT<50&&ZEN_PNT>40"``

is used to select observations with zenith angle range between 40 and 50 degrees.

The parameter ``datapath`` is only queried if the environment variable $VHEFITS
is not set.

The output of this script is an ASCII file containing observation IDs that
fulfil criteria specified by the user. This file can be used as input for
:doc:`csiactobs`. The script operates on an observation index file that needs
to be provided along with IACT data. `See here <http://gamma-astro-data-formats.readthedocs.org/en/latest/index.html>`__
how IACT data should be structured. To investigate additional parameters which
are available for selection, have a look `here <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/obs_index/index.html>`__.


General parameters
------------------

``datapath [string]``
    Path were data are located.

``prodname [string]``
    Name of FITS production (Run :doc:`csiactdata` to view your options).

``outfile [file]``
    Output runlist ASCII file.

``ra [real]``
    Right ascension of search region (invoke "INDEF" to avoid spatial selection).

``dec [real]``
    Declination of search region (invoke "INDEF" to avoid spatial selection).
    
``rad [real]``
    Radius of search region (deg) (invoke "INDEF" or negative values to avoid
    spatial selection)

``(min_qual = 0) [integer]``
    Minimum data quality (0=best, 1=medium, 2=bad).

``(expression = NONE) [string]``
    Additional expression.

``(master_index = master.json) [string]``
    Name of master index file.


Standard parameters
-------------------

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output runlist should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csfindobs.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`csiactobs`