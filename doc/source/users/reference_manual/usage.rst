Basic ctools parameter usage
============================

Automatic parameters
--------------------

``par_name [ = value ] <valid values> type``
    Where par_name is the name of the parameter, value is the default value,
    valid values are valid parameter values, separated by | characters, and
    type is the type of the parameter.
    The type is enclosed in square brackets.


Hidden parameters
-----------------

``(par_name = value) type``
    Where par_name is the name of the parameter, value is the default value,
    and type is the type of the parameter.
    The type is enclosed in square brackets.


Examples
--------

``infile [file]``
    Describes an automatic (queried) file-type parameter with no default value.
 	 	 
``irf = South_50h <North_50h|South_50h> [string]``
    Describes an automatic (queried) string-type parameter with two value
    options.
 	 	 
``(plot = yes) [boolean]``
    Describes a hidden bool-type parameter named plot, whose default value
    is yes (true).
 	 	 
.. note::

  To set hidden parameters at runtime, specify them explicitly at 
  the command line; for example:

  ``ctobssim start_time=900``
