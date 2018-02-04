.. _sec_introduction:

Introduction
============

ctools is a software package developed for the scientific analysis of 
Cherenkov Telescope Array (CTA) data or any other Imaging Air Cherenkov 
Telescope.

In this User Manual you will now learn how to use the ctools package.
If you encounter any problems, please read the :ref:`sec_help` section.
And you should read the :ref:`sec_develop` section if you would like to
contribute to the ctools developments.

Running ctools
--------------

ctools comprises a set of ftools-like binary executables, 
written in C++, with a command-line interface allowing for interactive
step-wise data analysis. To start a ctool from the command line, simply 
type the name of the ctool:

.. code-block:: bash

   $ ctobssim
   RA of pointing (degrees) (0-360) [83.63]

The ctool will query all user parameters and execute. User parameters can
also be specified as arguments separated by a whitespace to the ctool:

.. code-block:: bash

   $ ctobssim ra=83.63 dec=22.01 rad=5.0 tmin=0 tmax=1800 emin=0.1 emax=100.0 caldb=prod2 irf=South_0.5h \
              inmodel=$CTOOLS/share/models/crab.xml outevents=events.fits

If you need help about the usage of a ctool or a cscript, type the name of 
the ctool or cscript followed by the ``--help`` option:

.. code-block:: bash

   $ ctobssim --help

ctools includes also a Python module that declares every ctool as a
Python class. The following example illustrates the usage:

.. code-block:: bash

   $ python
   >>> import ctools
   >>> sim = ctools.ctobssim()
   >>> sim["ra"]  = 83.63
   >>> sim["dec"] = 22.01
   >>> sim.execute()
   Radius of FOV (degrees) (0-180) [5.0]

ctools comprise also cscripts, which are Python scripts that have the same 
user interface as the ctools binaries. cscripts can be started from the 
command line or used as Python classes, as illustrated below:

.. code-block:: bash

   $ python
   >>> import cscripts
   >>> pull = cscripts.cspull()

User parameters and parameter files
-----------------------------------

Each ctool and cscript has a defined list of user parameters. You can
find a complete description `here <../reference_manual/reference.html>`_.

User parameters are stored in par files. A ctools installation comprises a set
of default par files stored at ``$CTOOLS/syspfiles``. When you run a tool/script
the latest parameter values will be stored in a copy of the file at
``$HOME/pfiles`` (e.g., ``ctlike`` will save there a file called ``ctlike.par``).

When a tool/script is executed again, it will look for a par file first in
``$HOME/pfiles``, and propose you to use the latest values stored there. If no
par file is found there, it will use the default in ``$CTOOLS/syspfiles``.

If you update your ctools installation you may get an error such as

.. code-block:: bash
		
   *** ERROR encounterted in the execution of ctlike. Run aborted ...
   *** ERROR in GApplicationPars::load(GFilename&,
   std::vector<std::string>&): Invalid command line parameter
   encountered (invalid parameter name "XX"): XX=2

Remove the stale file ``ctlike.par`` from ``$HOME/pfiles`` and run again to
solve the issue.

