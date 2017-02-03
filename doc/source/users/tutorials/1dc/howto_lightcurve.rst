.. _1dc_howto_ligthcurve:

How to generate a light curve?
------------------------------

You can generate a light curve using the :ref:`cslightcrv` script which
takes on input an
:ref:`event list <glossary_eventlist>`
or a
:ref:`observation definition XML file <glossary_obsdef>`.
:ref:`cslightcrv` divides the event list into a number of intervals and
performs a maximum likelihood analyse for each of the intervals.
The results are writted into a light curve FITS file.

The following example illustrates how you can generate a light curve for
``Src001`` for the time interval from 51552.5 to 51552.6 (MJD). Ten time
intervals were chosen in this example.

.. code-block:: bash

   $ cslightcrv
   Input event list or observation definition XML file [events.fits] obs.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models.xml
   Source name [Crab] Src001
   Algorithm for defining time bins (FILE|LIN|GTI) [GTI] LIN
   Lightcurve start time (MJD) [51544.5] 51552.5
   Lightcurve stop time (MJD) [51544.6] 51552.6
   Number of time bins (1-10000) [5] 10
   Number of energy bins for binned (0=unbinned) (0-100) [0]
   Lower energy limit of events (TeV) [0.1]
   Upper energy limit of events (TeV) [100.0]
   Output light curve file [lightcurve.fits]

The resulting light curve is shown in the figure below.

.. figure:: howto_lightcurve.png
   :width: 600px
   :align: center

   *Light curve from 51552.5 to 51552.6 (MJD) for Src001*

.. note::
   The plot was created using the ``show_lightcurve.py`` script that is
   located in the ctools example folder. The example script requires the
   `matplotlib <http://matplotlib.org>`_ Python module for display.
   You may reproduce the plot by typing

   .. code-block:: bash

      $ $CTOOLS/share/examples/python/show_lightcurve.py lightcurve.fits
