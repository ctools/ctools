.. _1dc_first_select_obs:

Selecting the relevant observations
-----------------------------------

  .. admonition:: What you will learn

     You will learn how to **select a subset of observations** from the Data
     Challenge data for your analysis.

     We recommended to analyse the data region by region since the analysis of
     the Data Challenge data in a single shot will be very very time consuming.

Let's assume that you want to analyse the central region of our Galaxy using
the data obtained during the Galactic Plane Survey.

As first step you should create a folder that will contain your analysis
results and any intermediate data products. In this example we name the
folder ``my_first_analysis``, and create and step into it using

.. code-block:: bash

   $ mkdir my_first_analysis
   $ cd my_first_analysis

For the following, make sure that you have set the ``CTADATA`` and ``CALDB``
environment variables as desribed :ref:`here <1dc_getting_data>`.

The first step of the analysis consists in selecting the observations from the
:ref:`observation definition file <glossary_obsdef>`
``obs_gps_baseline.xml`` that have pointing directions close to the Galactic
Centre.
You do this with the :ref:`csobsselect` script by typing:

.. code-block:: bash

   $ csobsselect
   Input event list or observation definition XML file [obs.xml] $CTADATA/obs/obs_gps_baseline.xml
   Pointing selection region shape (CIRCLE|BOX) [CIRCLE]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] GAL
   Galactic longitude of selection centre (deg) (0-360) [184.56] 0.0
   Galactic latitude of selection centre (deg) (-90-90) [-5.79] 0.0
   Radius of selection circle (deg) (0-180) [5.0] 3.0
   Output observation definition XML file [outobs.xml] obs.xml

You selected ``CIRCLE`` as the shape of the pointing selection region,
specified the Galactic Centre as the centre of the selection region, and
selected all observations with pointing directions within 3 degrees from
the Galactic Centre from the ``obs_gps_baseline.xml`` file.
The selected observations are written into the file ``obs.xml`` that will be
created in the ``my_first_analysis`` folder.

You may use the :ref:`csobsinfo` script to display a summary of the selected
observations into the console:

.. code-block:: bash

   $ csobsinfo debug=yes
   Input event list, counts cube, or observation definition XML file [obs.xml]
   Output DS9 region file [ds9.reg]

.. note::

   The ``debug=yes`` attribute instructs :ref:`csobsinfo` to direct the log
   file output also to the console. Duplication of log file output into
   the console using the ``debug=yes`` attribute works for any tool or script.

The output of :ref:`csobsinfo` is shown below. :ref:`csobsselect` selected
57 observations from the Galactic Plane Survey data containing a total of
6 207 354 events.

.. code-block:: bash

  2017-06-02T08:10:00: +=========+
  2017-06-02T08:10:00: | Summary |
  2017-06-02T08:10:00: +=========+
  2017-06-02T08:10:00: === Observations ===
  2017-06-02T08:10:00:  Unbinned observations .....: 57
  2017-06-02T08:10:00:  Binned observations .......: 0
  2017-06-02T08:10:00: === Events ===
  2017-06-02T08:10:00:  Number of events ..........: 6207354
  2017-06-02T08:10:00:  Number of bins ............: 0
  2017-06-02T08:10:00: === Pointings ===
  2017-06-02T08:10:00:  Mean offset angle .........: Unknown
  2017-06-02T08:10:00:  Mean zenith angle .........: 0.00 deg
  2017-06-02T08:10:00:  Mean azimuth angle ........: 0.00 deg
  2017-06-02T08:10:00: === Energy range ===
  2017-06-02T08:10:00:  Minimum energy ............: 30 GeV
  2017-06-02T08:10:00:  Maximum energy ............: 160 TeV
  2017-06-02T08:10:00: === Time range ===
  2017-06-02T08:10:00:  MJD (days) ................: 59235.500 - 59276.921
  2017-06-02T08:10:00:  UTC .......................: 2021-01-21T11:58:51 - 2021-03-03T22:04:51
  2017-06-02T08:10:00:  MET (seconds) .............: 664502400.000 - 668081160.000
  2017-06-02T08:10:00:  Total ontime ..............: 102600.00 s = 1710.00 min = 28.50 h
  2017-06-02T08:10:00:  Total livetime ............: 100548.00 s = 1675.80 min = 27.93 h

The resulting
:ref:`observation definition file <glossary_obsdef>`
will look as follows:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="GPS" id="120380" instrument="CTA">
       <parameter name="EventList" file="/Users/jurgen/analysis/cta/dc/1dc/1dc.final/validation/1dc.south/data/baseline/gps/gps_baseline_120380.fits" />
       <parameter name="Calibration" database="1dc" response="South_z20_50h" />
     </observation>
     <observation name="GPS" id="120381" instrument="CTA">
       <parameter name="EventList" file="/Users/jurgen/analysis/cta/dc/1dc/1dc.final/validation/1dc.south/data/baseline/gps/gps_baseline_120381.fits" />
       <parameter name="Calibration" database="1dc" response="South_z20_50h" />
     </observation>
     ...
     <observation name="GPS" id="121177" instrument="CTA">
       <parameter name="EventList" file="/Users/jurgen/analysis/cta/dc/1dc/1dc.final/validation/1dc.south/data/baseline/gps/gps_baseline_121177.fits" />
       <parameter name="Calibration" database="1dc" response="South_z20_50h" />
     </observation>
   </observation_list>

Each ``<observation>`` element corresponds to one observation that is identified
by a ``name`` attribute and a unique identifier attribute.
An ``<observation>`` element contains two parameters:
the ``EventList`` parameter that specifies the name of the corresponding event
file, and
the ``Calibration`` parameter that specifies the
:ref:`instrument response function <glossary_irf>` that applies to the
event file.
