.. _start_selecting:

Selecting event data
--------------------

  .. admonition:: What you will learn

     You will learn how to use the :ref:`ctselect` tool to **select relevant
     events from the simulated event data**.

As next step you may want to select some events from the simulated event data.
You may want to restrict the events to a region of interest (ROI) of 3 degrees
around the pointing direction, you may want to extract a 30 minutes time slice
from the data, and you may want to limit the energy range of the events
to 0.1 - 100 TeV.
You do this event selection using the :ref:`ctselect` tool as follows:

.. code-block:: bash

   $ ctselect
   Input event list or observation definition XML file [events.fits]
   Radius of ROI around pointing or specified RA/DEC (0-180) [3.0]
   Start time (UTC string, JD, MJD or MET in seconds) [NONE] 2020-01-01T00:10:00
   Stop time (UTC string, JD, MJD or MET in seconds) [NONE] 2020-01-01T00:40:00
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Output event list or observation definition XML file [selected_events.fits]

This created a new event list FITS file ``selected_events.fits``
that contains the selected data.
Below an excerpt of the ``ctselect.log`` log file that illustrates the
event selection:

.. code-block:: none

   2018-01-24T14:07:07: +===================+
   2018-01-24T14:07:07: | Input observation |
   2018-01-24T14:07:07: +===================+
   2018-01-24T14:07:07: === GObservations ===
   2018-01-24T14:07:07:  Number of observations ....: 1
   2018-01-24T14:07:07:  Number of models ..........: 0
   2018-01-24T14:07:07:  Number of observed events .: 202210
   2018-01-24T14:07:07:  Number of predicted events : 0
   2018-01-24T14:07:07:
   2018-01-24T14:07:07: +=================+
   2018-01-24T14:07:07: | Event selection |
   2018-01-24T14:07:07: +=================+
   2018-01-24T14:07:07: === CTA observation (id=000001) ===
   2018-01-24T14:07:07:  Input filename ............: events.fits
   2018-01-24T14:07:07:  Event extension name ......: EVENTS
   2018-01-24T14:07:07:  GTI extension name ........: GTI
   2018-01-24T14:07:07: === Event selection ===
   2018-01-24T14:07:07:  Time range (MJD) ..........: 58849.0077451852 - 58849.0285785185 days
   2018-01-24T14:07:07:  Time range (UTC) ..........: 2020-01-01T00:10:00 - 2020-01-01T00:40:00
   2018-01-24T14:07:07:  Time range (MET) ..........: 631109469.184 - 631111269.184 seconds
   2018-01-24T14:07:07:  Selected energy range .....: 0.1 - 100 TeV
   2018-01-24T14:07:07:  Requested RoI .............: Centre(RA,DEC)=(83.63, 22.51) deg, Radius=3 deg
   2018-01-24T14:07:07:  RoI of data ...............: Centre(RA,DEC)=(83.63, 22.51) deg, Radius=5 deg
   2018-01-24T14:07:07:  Selected RoI ..............: Centre(RA,DEC)=(83.63, 22.51) deg, Radius=3 deg
   2018-01-24T14:07:07:  cfitsio selection .........: TIME >= 631109469.18399978 && TIME <= 631111269.18400002 && ENERGY >= 0.10000000 && ENERGY <= 100.00000000 && ANGSEP(83.630000,22.510000,RA,DEC) <= 3.000000
   2018-01-24T14:07:07:  FITS filename .............: /var/tmp/tmp.0.aFKvXb[EVENTS][TIME >= 631109469.18399978 && TIME <= 631111269.18400002 && ENERGY >= 0.10000000 && ENERGY <= 100.00000000 && ANGSEP(83.630000,22.510000,RA,DEC) <= 3.000000]
   2018-01-24T14:07:07:
   2018-01-24T14:07:07: +====================+
   2018-01-24T14:07:07: | Output observation |
   2018-01-24T14:07:07: +====================+
   2018-01-24T14:07:07: === GObservations ===
   2018-01-24T14:07:07:  Number of observations ....: 1
   2018-01-24T14:07:07:  Number of models ..........: 0
   2018-01-24T14:07:07:  Number of observed events .: 21974
   2018-01-24T14:07:07:  Number of predicted events : 0

From the 202210 events in the input event list 21974 events were selected.

.. note::
   Event selection by :ref:`ctselect` in a given dimension can be avoided
   by specifying ``NONE`` (or ``INDEF``, ``UNDEF`` or ``UNDEFINED``) for
   one of the selection parameters. This is recommended rather than
   specifying the same selection previosuly applied to the data. The following run for example only
   performs a selection by time but no selection by ROI or energy:

   .. code-block:: bash

      $ ctselect
      Input event list or observation definition XML file [events.fits]
      Radius of ROI around pointing or specified RA/DEC (degrees) (0-180.) [3.] NONE
      Start time (UTC string, JD, MJD or time in seconds) [2020-01-01T00:10:00]
      Stop time (UTC string, JD, MJD or time in seconds) [2020-01-01T00:40:00]
      Lower energy limit (TeV) [0.1] NONE
      Output event list or observation definition XML file [selected_events.fits]
