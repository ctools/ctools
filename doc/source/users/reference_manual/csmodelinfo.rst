.. _csmodelinfo:

csmodelinfo
===========

Dumps information about a model container into log file or on screen.


Synopsis
--------

This script provides detailled information about a given model container.
It presents the number of parameters, the number of free parameters (sub-
divided by spectral, spatial and temporal type). The tool further creates
information about any parameter that is stuck at the parameter boundaries.
This might be quite useful for debugging a fit. Finally the tool also list
the TS values of each source that was fitted, in order have a quick overview
on these quantities on the screen.

In addition the tool is also capable of exporting the input model XML container
to a region file which is compatible to be read with ds9 or other displying
tools. To create a region file, the hidden parameter "ds9file" needs to be
specified. There are quite a few further hidden parameters to steer e.g.
color, font, etc of the region file output.


General parameters
------------------

``inmodel [file]``
    Event list, counts cube, or observation definition file

``outds9file [file]``
    Output DS9 region file containing soucre positions

``(pnt_type=cross) <circle|box|diamond|cross|x|arrow|boxcircle> [string]``
    Marker type for point sources

``pnt_mark_size [integer]``
    Marker size for point sources

``(show_labels=yes) [boolean]``
    Add source labels?
    
``(width=2) [integer]``
    Line width for regions
    
 ``(font=helvetica) <helvetica|times|courier> [string]``
    Font for source labels

``(fontsize=12) [integer]``
    Font size for source labels

``(fonttype=normal) <normal|bold> [string]``
    Use normal or bold font?

``(fonttype2=roman) <roman|italic> [string]``
    Use roman or italic font?

``(show_ext_type=yes) [boolean]``
    Show type of extended model in source name?

``(free_color=green) [string]``
    Color for sources with free parameters (any ds9 color or hex code)

``(fixed_color=magenta) [string]``
    Color for sources without free parameters (any ds9 color or hex code)
    
    
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
    Specifies whether existing output files should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csmodelinfo.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

None
