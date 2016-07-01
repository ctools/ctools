.. _csmodelmerge:

csmodelmerge
============

Merges an arbitrary number of model XML files into a single file.


Synopsis
--------

This script merges several model XML file into a single model XML file. The
files are passed via the ``inmodels`` parameter. There are 4 options to specify
a list of files that should be merged:

* a space-separated list of file names (e.g. ``model1.xml model2.xml``)
* a semi-colon separated list of file names (e.g. ``model1.xml;model2.xml``)
* a wildcard string (e.g. ``mymodel*.xml``)
* an ASCII file containing the file names, one file per line (e.g. ``@mymodels.txt``)

Note that the "@" needs to be specified in case of an ASCII file name. It is
also important to note that models with the same name cannot be appended. The
tool will throw and exception if the provided file contain models with the same
name.


General parameters
------------------

``inmodels [string]``
    Input model XML files

``outmodel [file]``
    Output model file
    
    
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
    Specifies whether an existing output model file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csmodelmerge.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

None