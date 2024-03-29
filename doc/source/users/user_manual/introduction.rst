.. _sec_introduction:

Introduction
============

ctools is a software package developed for the scientific analysis of
gamma-ray astronomy data with emphasis on data from the Cherenkov Telescope
Array (CTA) data or any other Imaging Air Cherenkov Telescope.

In this User Manual you will now learn how to use the ctools package.
If you encounter any problems, please read the :ref:`sec_help` section.
And you should read the :ref:`sec_develop` section if you would like to
contribute to the ctools developments.

The ctools software package is based on `GammaLib <http://cta.irap.omp.eu/gammalib/index.html>`_,
a C++ API and Python module that provides support for the scientific analysis of
astronomical gamma-ray data. For more information on
`GammaLib <http://cta.irap.omp.eu/gammalib/index.html>`_
please read the `GammaLib User documentation <http://cta.irap.omp.eu/gammalib/users/index.html>`_.


Running ctools
--------------

ctools comprises a set of ftools-like binary executables written in C++
and scripts written in Python. Binary executables start with ``ct`` (for
example :ref:`ctobssim` for the simulation of data), Python scripts start
with ``cs`` (for example :ref:`cslightcrv` for the generation of a light
curve). While a binary executable is called a *ctool*, a Python script is
called a *cscript*. From a user perspective, ctools and cscripts behave
identically. Both have a command-line interface allowing for interactive
step-wise analysis of the data. From now on we will simply say *tool* when
considering either a ctool or a cscript.

To use a tool from the command line, simply type the name of the tool. For
example, to simulate CTA data you can type:

.. code-block:: bash

   $ ctobssim
   RA of pointing (degrees) (0-360) [83.63]

The tool will query all user parameters and execute. User parameters can
also be specified as arguments separated by a whitespace to the tool:

.. code-block:: bash

   $ ctobssim ra=83.63 dec=22.01 rad=5.0 tmin=0 tmax=1800 emin=0.1 emax=100.0 caldb=prod2 irf=South_0.5h \
              inmodel=$CTOOLS/share/models/crab.xml outevents=events.fits

If you need help about the usage of a tool, type the name of the tool followed
by the ``--help`` option:

.. code-block:: bash

   $ ctobssim --help

ctools includes also Python modules that declare every tool as a Python class.
ctools are available through the ``ctools`` module, cscripts through the
``cscripts`` module. The following example illustrates how data can be simulated
from Python:

.. code-block:: python

   >>> import ctools
   >>> sim = ctools.ctobssim()
   >>> sim["ra"]  = 83.63
   >>> sim["dec"] = 22.01
   >>> sim.execute()
   Radius of FOV (degrees) (0-180) [5.0]

And here an example for generating a light curve from Python:

.. code-block:: python

   >>> import cscripts
   >>> lightcrv = cscripts.cslightcrv()
   >>> lightcrv.execute()


User parameters and parameter files
-----------------------------------

Each ctool and cscript has a defined list of user parameters. You can
find a complete description of the parameters of all tools
`here <../reference_manual/index.html>`_.

User parameters are stored in par files. A ctools installation comprises a set
of default par files stored at ``$CTOOLS/syspfiles``. When you run a tool/script
the latest parameter values will be stored in a copy of the file at
``$HOME/pfiles`` (e.g., ``ctlike`` will save there a file called ``ctlike.par``).

When a tool/script is executed again, it will look for a par file first in
``$HOME/pfiles``, and propose you to use the latest values stored there. If no
par file is found there, it will use the default in ``$CTOOLS/syspfiles``.

If you delete files in ``$HOME/pfiles`` the latest values will be lost, and upon
execution of a tool a new copy of the par file will be stored in the folder. If
for some reason a parameter file got corrupt, simply delete it from ``$HOME/pfiles``
and start again.

There are two types of user parameters: those that will be queried when starting
a tool and those that are hidden and not queried. Hidden parameters serve to
define default values that normally need not to be changed, but exposing the
parameters in the interface allows to change them by the user for fine tuning
of a tool. An example for a hidden parameter of :ref:`ctobssim` is the ``seed``
parameter to set the initial seed value of the random number generator. To set
a hidden parameter its value has to be specified on the command line. For
example

.. code-block:: bash

   $ ctobssim seed=41

will run :ref:`ctobssim` with a seed value of 41. Multiple hidden parameters
specified on the command line need to be separated by a white space. In the
Python interface, values of hidden parameters are specified in the same way as
values of queried parameters, e.g.

.. code-block:: python

   >>> import ctools
   >>> sim = ctools.ctobssim()
   >>> sim["seed"] = 41
   ...


Log files
---------

Each ctool and cscript will produce a log file that is written in the working
directory. The log file is in ASCII file format and can be read by any editor.
By default, the name of the ASCII file is the name of the tool or script,
suffixed by ``.log``. You can change the log file name using the hidden
``logfile`` parameter, for example

.. code-block:: bash

   $ ctobssim logfile=my_first_simulation.log

When calling a ctool or cscript from Python, no log file will be written by
default. The reason for this is that Python scripts are often used to build
ctools analysis pipelines and workflows, and one generally does not want that
such a script pollutes the workspace with log files. You can however instruct
a ctool or cscript to generate a log file by invoking the ``logFileOpen()``
method before running the tool, for example

.. code-block:: python

   >>> import ctools
   >>> sim = ctools.ctobssim()
   ...
   >>> sim.logFileOpen()
   >>> sim.run()


Additional modules
------------------

From version 2.0.0 on ctools includes an additional module comprising scripts
for COMPTEL/CGRO data analysis. All COMPTEL scripts start with ``com`` (for example
:ref:`comlixfit` for the fitting of data using the SRCLIX algorithm) and are used
in the same way as the regular ctools or cscripts. All COMPTEL scripts are also
available from Python through a dedicated ``comscripts`` module:

.. code-block:: python

   >>> import comscripts
   >>> fit = comscripts.comlixfit()

