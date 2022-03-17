.. _tutorials_fermi:

Fermi/LAT
=========

  .. admonition:: What you will learn

     You will learn how to analyse data from the LAT telescope aboard NASA's
     `Fermi <https://fermi.gsfc.nasa.gov/ssc/>`_ satellite with ctools.

You can use ctools to analyse data from the LAT telescope aboard NASA's
`Fermi <https://fermi.gsfc.nasa.gov/ssc/>`_ satellite.
To do this you need to be familar with the Fermi-LAT Science Tools that you
will use to prepare the Fermi-LAT data for the analysis. You can download the
Fermi-LAT Science Tools from the
`Fermi Science Support Center (FSSC) <https://fermi.gsfc.nasa.gov/ssc/>`_
which provides also more detailed information about how to prepare the data.
What you need for a given observation are the following files:

* ``srcmaps.fits`` - source maps file that include a counts cube and maps of the
  diffuse model components
* ``ltcube.fits`` - livetime cube file that corresponds to the Good Time Intervals
  of the data
* ``expmap.fits`` - exposure map file

To demonstrate how the Fermi-LAT analysis works with ctools, here is a short
tutorial that is based on a simple analysis of the Vela pulsar.

.. toctree::
   :maxdepth: 1

   fermi_prepare
   fermi_fitting
   fermi_butterfly
   fermi_spectrum
   fermi_tsmap
   fermi_ulimit
   fermi_errors
