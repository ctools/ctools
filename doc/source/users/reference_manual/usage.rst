.. _sec_usage:

ctools parameter interface
==========================

For a generic introduction to ctools and user parameters see
`here <../user_manual/introduction.html>`__.

Automatic parameters
--------------------

``par_name <first valid value|second valid value|...> [type]``
    Where ``par_name`` is the name of the parameter, ``first valid value`` and
    ``second valid value`` are valid parameter values, separated by ``|`` characters,
    and ``type`` is the type of the parameter.


Hidden parameters
-----------------

``(par_name = value) [type]``
    Where ``par_name`` is the name of the parameter, ``value`` is the default value,
    and ``type`` is the type of the parameter.


Examples
--------

``infile [file]``
    Describes an automatic (queried) file-type parameter.

``irf <North_50h|South_50h> [string]``
    Describes an automatic (queried) string-type parameter with two value
    options.

``(plot = yes) [boolean]``
    Describes a hidden bool-type parameter named plot, whose default value
    is yes (true).

.. note::

   You can set hidden parameters at runtime by specifying them explicitly on 
   the command line; for example:

   ``ctobssim start_time=900``
