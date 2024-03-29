.. _um_iact_data:

IACT data
---------

For Imaging Air Cherenkov Telescopes (IACTs), such as CTA, H.E.S.S., MAGIC or
VERITAS, ctools operates on event lists. An event list is a table where each row
corresponds to a registered event and each column corresponds to a property
of the event. Required properties are the event identifier, the reconstructed
event direction, the reconstructed event energy and the trigger time. Event
lists are provided as FITS binary tables. An example for a minimal event list
is shown below.

.. figure:: iact_eventlist.png
   :width: 70%
   :align: center

Metadata, such as the telescope pointing direction, the live time, or any other
information that may be relevant for data processing, is stored in the FITS
header of the event list binary table. An example for the metadata that is
included in the header of the H.E.S.S. Collaboration event lists is shown
below.

.. figure:: iact_eventheader.png
   :width: 70%
   :align: center

Every event list needs to be accompanied by a table with so-called
:ref:`Good Time Intervals (GTIs) <glossary_gti>` where each row corresponds
to a time interval of continuous data taking, and the two columns correspond
to the start and stop times of these time intervals. An example for the GTIs
of the H.E.S.S. event lists comprising a single GTI is shown below.

.. figure:: iact_gti.png
   :width: 40%
   :align: center

Note that times are given in so-called Mission-Elapsed Time (MET) in seconds,
and that the zero-point of the MET is specified by the ``MJDREFI`` and
``MJDREFF`` keywords as the integer and fractional part of a Modified Julian Date
(MJD) in the header of the GTI. The same holds also for the trigger times in
the event list.

For further information on the data format we recommend to read the
`Data formats for gamma-ray astronomy <https://gamma-astro-data-formats.readthedocs.io/en/latest/>`_
document.
