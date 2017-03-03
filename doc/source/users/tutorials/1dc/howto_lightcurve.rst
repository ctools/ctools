.. _1dc_howto_ligthcurve:

How to generate a light curve?
------------------------------

You can generate a light curve using the :ref:`cslightcrv` script which
takes on input an
:ref:`event list <glossary_eventlist>`
or an
:ref:`observation definition file <glossary_obsdef>`.
Light curves cannot be generated from counts cubes since the time information
is dropped in the generation of the counts cube.

:ref:`cslightcrv` divides the event list into a number of intervals and
performs a maximum likelihood analyse for each of the intervals.
The analysis results are writted into a light curve FITS file.

The following example illustrates how you can generate a light curve for
``Src001`` for the time interval from 59236.810 to 59245.742 (MJD),
corresponding to the dates from 2021-01-22 (19:26:51) to 2021-01-31 (17:48:51)
(to convert dates to MJD and vice versa we suggest to use the
`xTime Date and Time Conversion Utility <https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime/xTime.pl>`_).
The time interval is linearly divided into 100 time bins.

.. code-block:: bash

   $ cslightcrv
   Input event list or observation definition XML file [events.fits] obs_selected.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models.xml
   Source name [Crab] Src001
   Algorithm for defining time bins (FILE|LIN|GTI) [GTI] LIN
   Lightcurve start time (MJD) [51544.5] 59236.810
   Lightcurve stop time (MJD) [51544.6] 59245.742
   Number of time bins (1-10000) [5] 100
   Number of energy bins for binned (0=unbinned) (0-100) [0]
   Lower energy limit of events (TeV) [0.1]
   Upper energy limit of events (TeV) [100.0]
   Output light curve file [lightcurve.fits]

The resulting light curve is shown in the figure below. The gaps in the light
curve correspond to time periods during which ``Src001`` has not been observed.

.. figure:: howto_lightcurve.png
   :width: 600px
   :align: center

   *Light curve from 59236.810 to 59245.742 (MJD) for Src001*

.. note::
   The plot was created using the ``show_lightcurve.py`` script that is
   located in the ctools example folder. The example script requires the
   `matplotlib <http://matplotlib.org>`_ Python module for display.
   You may reproduce the plot by typing

   .. code-block:: bash

      $ $CTOOLS/share/examples/python/show_lightcurve.py lightcurve.fits
