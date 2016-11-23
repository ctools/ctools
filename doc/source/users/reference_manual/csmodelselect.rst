.. _csmodelselect:

csmodelselect
=============

Selects models from a model definition XML file.


Synopsis
--------

This script selects all models from a model definiton XML file that spatially
overlap with the Region of Interests (RoIs) of a list of observations and writes
these models into a new model definition XML file. In addition to the spatial
selection, the tool allows selecting only models with a photon flux and/or a
Test Statistic value above a given limit. The script also allows to specify
whether the source positions and/or source shapes should be fit.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file

``inmodel [file]``
    Input model definition XML file

``outmodel [file]``
    Output model definition XML file

``(roilimit = 4.5) [real]``
    Maximum RoI radius (degrees)

``(roimargin = 0.1) [real]``
    Radial margin to be added to RoIs (degrees)

``(ethres = 0.1) [real]``
    Energy threshold for source flux selection (TeV)

``(fluxlimit = 1.0e-12) [real]``
    Minimum source flux for selection (ph/cm2/s)

``(tslimit = 10.0) [real]``
    Minimum Test Statistic for selection

``(fit_pos = yes) [boolean]``
    Fit source positions?

``(fit_shape = no) [boolean]``
    Fit source shapes?

    
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
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csmodelselect.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`csobsselect`
