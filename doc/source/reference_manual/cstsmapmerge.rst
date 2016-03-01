.. _cstsmapmerge:

cstsmapmerge
============

Merges slices of splitted tsmap computations.


Synopsis
--------

In case a TS map has been computed using the option to split the calculations
into several jobs, this tool is required to merge the output into one final
map again. Running :doc:`cttsmap` with the hidden parameters ``binmin`` and
``binmax`` specified will result in several sliced TS maps. This script merges
them into one FITS file. The input slices can be passed via the ``inmaps``
parameter. There are 4 options to specify a list of files that should be merged:

* a space-separated list of file names (e.g. tsmap1.fits tsmap2.fits)
* a comma-separated list of file names (e.g. tsmap1.fits,tsmap2.fits)
* a wildcard string (e.g. tsmap*.xml)
* an ASCII file containing the file names, one file per line (e.g. @mymaps.txt)

Note that the "@" needs to be specified in case of an ASCII file name. 

General parameters
------------------

``inmaps [string]``
    Input TS map FITS files

``outmap [file]``
    Output TS map FITS file
    
``(overwrite=yes) [boolean]``
    Overwrite previously filled bins?

``(delete=no) [boolean]``
    Delete merged files?
    
    
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
    Specifies whether an existing output TS map file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = cstsmapmerge.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:ref:`cttsmap`
