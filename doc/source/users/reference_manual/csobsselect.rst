.. _csobsselect:

csobsselect
===========

Selects observations from an observation definition XML file.


Synopsis
--------

This script selects all observations from an observation definiton XML file that
have pointings within a specified selection region and writes these observations
into a new observation definition XML file. Possible selection regions are a
circle and a box, either in Galactic or Celestial coordinates.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file

``outobs [file]``
    Output observation definition XML file

``pntselect <CIRCLE|BOX> [string]``
    Pointing selection region shape

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic)

``ra [real]``
    Right Ascension of pointing selection centre (J2000, in degrees)
 	 	 
``dec [real]``
    Declination of pointing selection centre (J2000, in degrees)

``glon [real]``
    Galactic longitude of pointing selection centre (degrees)
 	 	 
``glat [real]``
    Galactic latitude of pointing selection centre (degrees)

``rad [real]``
    Radius of selection circle (degrees)

``width [real]``
    Width of selection box (degrees)

``height [real]``
    Height of selection box (degrees)

    
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
    Specifies whether an existing output model file should be overwritten
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn")

``(logfile = csobsselect.log) [filename]``
    Log filename


Related tools or scripts
------------------------

:doc:`csmodelselect`

