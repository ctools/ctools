.. _comptel_selection:

Select data from database
-------------------------

  .. admonition:: What you will learn

     You will learn how to select COMPTEL data for an analysis with ctools.


The first step of a COMPTEL science analysis consists in the selection
of the relevant viewing periods. You do this with the :ref:`comobsselect`
script. In the following example, all viewing periods with pointing directions
within 30 degrees of the Crab and with dates between April and May 1991
are selected.

.. code-block:: bash

   $ comobsselect
   Pointing selection region shape (CIRCLE|BOX) [CIRCLE]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [GAL]
   Galactic longitude of selection centre (deg) (0-360) [184.56]
   Galactic latitude of selection centre (deg) (-90-90) [-5.79]
   Radius of selection circle (deg) (0-180) [30]
   Start time (UTC string, JD, MJD or MET in seconds) [NONE] 1991-05-01T00:00:00
   Stop time (UTC string, JD, MJD or MET in seconds) [NONE] 1991-06-01T00:00:00
   Output observation definition file [obs.xml]

Giving the limits on the start and stop time, only viewing period ``0001``
will be selected by :ref:`comobsselect`. The resulting
:ref:`observation definition file <glossary_obsdef>`
``obs.xml`` will have the following content:

.. code-block:: bash

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="CRAB" id="vp0001_0" instrument="COM">
       <parameter name="EVP" file="$COMDATA/phase01/vp0001_0/m28511_evp.fits" />
       <parameter name="TIM" file="$COMDATA/dbase/tim/vp0001_0_tim.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20037_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20039_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20041_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20045_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20048_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20050_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20054_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20058_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20064_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20066_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20068_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20071_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20073_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20078_oad.fits" />
       <parameter name="OAD" file="$COMDATA/phase01/vp0001_0/m20081_oad.fits" />
     </observation>
   </observation_list>

The file contains a single observation, composed of an event file (``EVP``),
a Good Time Interval file (``TIM``) and 15 Orbit Aspect data files (``OAD``).

.. note::

   If no limit on the observation dates should be applied then specify ``NONE``
   as the start or stop time.
