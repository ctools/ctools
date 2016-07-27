.. _sec_introduction:

Introduction
============

ctools is a software package developed for the scientific analysis of 
Cherenkov Telescope Array (CTA) data or any other Imaging Air Cherenkov 
Telescope. ctools comprises a set of ftools-like binary executables, 
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

In this User Manual you will now learn how to install and use the ctools
package.
If you encounter any problems, please read the :ref:`faq` and :ref:`issues`
sections. If you're desperate, you may need :ref:`help`.
And you should read the :ref:`develop` section if you would like to 
contribute to the ctools developments.
