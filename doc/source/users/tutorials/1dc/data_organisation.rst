.. _1dc_data_organisation:

Data organisation
-----------------

Layout
^^^^^^

The content of the ``1dc`` folder should be as follows:

.. code-block:: bash

   caldb/
   caldb/data
   caldb/data/cta
   caldb/data/cta/prod3b
   caldb/data/cta/prod3b/caldb.indx
   caldb/data/cta/prod3b/bcf
   ...
   data/
   data/gc_baseline_000001.fits
   data/gc_baseline_000002.fits
   ...
   data/gc_threshold_000001.fits
   data/gc_threshold_000002.fits
   ...
   data/
   data/gps_baseline_000001.fits
   data/gps_baseline_000002.fits
   ...
   data/gps_threshold_000001.fits
   data/gps_threshold_000002.fits
   ...
   data/egal_baseline_000001.fits
   data/egal_baseline_000002.fits
   ...
   data/egal_threshold_000001.fits
   data/egal_threshold_000002.fits
   ...
   data/obs_gc_baseline.xml
   data/obs_gc_threshold.xml
   data/obs_gps_baseline.xml
   data/obs_gps_threshold.xml
   data/obs_egal_baseline.xml
   data/obs_egal_threshold.xml
   data/obs_all_baseline.xml
   data/obs_all_threshold.xml
   models/
   models/models_gc.xml
   models/models_gps.xml
   models/models_egal.xml
   models/map_RXJ1713.fits.gz
   models/map_VelaJunior.fits.gz
   models/map_ics.fits.gz
   models/map_pi0.fits.gz


Instrument Response Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``caldb`` folder contains the
:ref:`Instrument Response Functions <glossary_irf>`
that are necessary for the analysis of the simulated CTA data.
The folder contains the ``prod3b`` response that should be used for the
:ref:`first CTA Data Challenge <glossary_1dc>`.
The following response functions are available:

 +-----------------------+-------+---------------+--------+----------+
 | Response name         | Site  | Configuration | Zenith | Duration |
 +=======================+=======+===============+========+==========+
 | ``South_z20_50h``     | South | Baseline      | 20 deg | 50 hours |
 +-----------------------+-------+---------------+--------+----------+
 | ``South_z40_50h``     | South | Baseline      | 40 deg | 50 hours |
 +-----------------------+-------+---------------+--------+----------+
 | ``South_z20_5h``      | South | Baseline      | 20 deg | 5 hours  |
 +-----------------------+-------+---------------+--------+----------+
 | ``South_z40_5h``      | South | Baseline      | 40 deg | 5 hours  |
 +-----------------------+-------+---------------+--------+----------+
 | ``South_z20_0.5h``    | South | Baseline      | 20 deg | 30 min   |
 +-----------------------+-------+---------------+--------+----------+
 | ``South_z40_0.5h``    | South | Baseline      | 40 deg | 30 min   |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_z20_50h``     | North | Baseline      | 20 deg | 50 hours |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_z40_50h``     | North | Baseline      | 40 deg | 50 hours |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_z20_5h``      | North | Baseline      | 20 deg | 5 hours  |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_z40_5h``      | North | Baseline      | 40 deg | 5 hours  |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_z20_0.5h``    | North | Baseline      | 20 deg | 30 min   |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_z40_0.5h``    | North | Baseline      | 40 deg | 30 min   |
 +-----------------------+-------+---------------+--------+----------+
 | ``South_TS_z20_50h``  | South | Threshold     | 20 deg | 50 hours |
 +-----------------------+-------+---------------+--------+----------+
 | ``South_TS_z40_50h``  | South | Threshold     | 40 deg | 50 hours |
 +-----------------------+-------+---------------+--------+----------+
 | ``South_TS_z20_5h``   | South | Threshold     | 20 deg | 5 hours  |
 +-----------------------+-------+---------------+--------+----------+
 | ``South_TS_z40_5h``   | South | Threshold     | 40 deg | 5 hours  |
 +-----------------------+-------+---------------+--------+----------+
 | ``South_TS_z20_0.5h`` | South | Threshold     | 20 deg | 30 min   |
 +-----------------------+-------+---------------+--------+----------+
 | ``South_TS_z40_0.5h`` | South | Threshold     | 40 deg | 30 min   |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_TS_z20_50h``  | North | Threshold     | 20 deg | 50 hours |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_TS_z40_50h``  | North | Threshold     | 40 deg | 50 hours |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_TS_z20_5h``   | North | Threshold     | 20 deg | 5 hours  |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_TS_z40_5h``   | North | Threshold     | 40 deg | 5 hours  |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_TS_z20_0.5h`` | North | Threshold     | 20 deg | 30 min   |
 +-----------------------+-------+---------------+--------+----------+
 | ``North_TS_z40_0.5h`` | North | Threshold     | 40 deg | 30 min   |
 +-----------------------+-------+---------------+--------+----------+

.. warning::
   The **50 hours**
   :ref:`Instrument Response Functions <glossary_irf>`
   were used for the **simulation** of the
   :ref:`first CTA Data Challenge <glossary_1dc>`
   data. Please use only these response functions for the analysis. If you use
   :ref:`Observation Definition Files <glossary_obsdef>`
   for the analysis (see below) the appropriate 50 hours response functions
   will be used automatically.


Event data
^^^^^^^^^^

The ``data`` folder contains the calibrated, reconstructed and background
reduced event data that were procuded for the
:ref:`first CTA Data Challenge <glossary_1dc>`
and that were stored into FITS files.
Each event file contains the events for an
:ref:`observation <glossary_obs>`
(or run) of 30 minutes duration and comprises an
:ref:`event list <glossary_eventlist>`
and a
:ref:`Good Time Intervals <glossary_gti>`
binary table extension (see figure below).

.. figure:: event_file.png
   :width: 600px
   :align: center

   *Structure of an event file*

The header of the ``EVENTS`` table contains information about the
:ref:`observation <glossary_obs>`
such as
the start and stop date and time,
the duration and livetime of the observation, and
the pointing direction in Right Ascension and Declination (see figure below).

.. figure:: event_header.png
   :width: 500px
   :align: center

   *Header of an event list*

.. note::
   The pointing direction during an observation is fixed.

.. warning::
   Only the following header keywords in the ``EVENTS`` table have meaningful
   values:

   * ``DATE_OBS`` - start date of observation (UTC)
   * ``TIME_OBS`` - start time of observation (UTC)
   * ``DATE_END`` - end date of observation (UTC)
   * ``TIME_END`` - end time of observation (UTC)
   * ``TSTART`` - start time of observation, counted from time reference (s)
   * ``TSTOP`` - stop time of observation, counted from time reference (s)
   * ``MJDREFI`` - integer part of time reference MJD (days)
   * ``MJDREFF`` - fractional part of time reference MJD (days)
   * ``TIMEUNIT`` - time unit
   * ``TIMESYS`` - time system
   * ``TIMEREF`` - time reference
   * ``TELAPSE`` - elapsed time (s)
   * ``ONTIME`` - exposure time (s)
   * ``LIVETIME`` - livetime (s)
   * ``DEADC`` - deadtime correction factor, livetime / exposure time
   * ``TIMEDEL`` - time resolution
   * ``RA_PNT`` - Right Ascension of pointing direction (deg)
   * ``DEC_PNT`` - Declination of pointing direction (deg)
   * ``RADECSYS`` - Coordinate system
   * ``EQUINOX`` - Coordinate epoch

   All remaining header keywords have arbitrary values and should not be
   used for the analysis.

Each row of the ``EVENTS`` table corresponds to a single event.
Each event is characterised by

 +--------------+-------------------------------------------+----------+
 | Column       | Meaning                                   | Unit     |
 +==============+===========================================+==========+
 | ``EVENT_ID`` | Event number in file                      | unitless |
 +--------------+-------------------------------------------+----------+
 | ``TIME``     | Time stamp, countered from time reference | s        |
 +--------------+-------------------------------------------+----------+
 | ``RA``       | Reconstructed Right Ascension             | deg      |
 +--------------+-------------------------------------------+----------+
 | ``DEC``      | Reconstructed Declination                 | deg      |
 +--------------+-------------------------------------------+----------+
 | ``ENERGY``   | Reconstructed energy                      | TeV      |
 +--------------+-------------------------------------------+----------+
 | ``DETX``     | Reconstructed camera X coordinate         | deg      |
 +--------------+-------------------------------------------+----------+
 | ``DETY``     | Reconstructed camera Y coordinate         | deg      |
 +--------------+-------------------------------------------+----------+

An example of an ``EVENTS`` table is shown below.

.. figure:: event_list.png
   :width: 600px
   :align: center

   *Content of an event list*

Event files exist for

 +-----------------+---------------+--------------------+------------------+
 | Survey          | Configuration | Prefix             | Energy range     |
 +=================+===============+====================+==================+
 | Galactic Centre | Baseline      | ``gc_baseline``    |                  |
 +-----------------+---------------+--------------------+------------------+
 | Galactic Plane  | Baseline      | ``gps_baseline``   | 30 GeV - 120 TeV |
 +-----------------+---------------+--------------------+------------------+
 | Extragalactic   | Baseline      | ``egal_baseline``  |                  |
 +-----------------+---------------+--------------------+------------------+
 | Galactic Centre | Threshold     | ``gc_threshold``   |                  |
 +-----------------+---------------+--------------------+------------------+
 | Galactic Plane  | Threshold     | ``gps_threshold``  |                  |
 +-----------------+---------------+--------------------+------------------+
 | Extragalactic   | Threshold     | ``egal_threshold`` |                  |
 +-----------------+---------------+--------------------+------------------+


Observation Definition Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The files ``obs_XXX_YYYYY.xml`` are so called
:ref:`Observation Definition Files <glossary_obsdef>`
that contain the information (or metadata) of a list of observations.
The files are plain ASCII files in XML format that can be inspected and
manipulated by any text editor.
The following files are included in the distribution:

 +----------------------------+--------------------------------------------------------+
 | File                       | Purpose                                                |
 +============================+========================================================+
 | ``obs_gc_baseline.xml``    | Galactic Centre Survey performed with baseline arrays  |
 +----------------------------+--------------------------------------------------------+
 | ``obs_gc_threshold.xml``   | Galactic Centre Survey performed with threshold arrays |
 +----------------------------+--------------------------------------------------------+
 | ``obs_gps_baseline.xml``   | Galactic Plane Survey performed with baseline arrays   |
 +----------------------------+--------------------------------------------------------+
 | ``obs_gps_threshold.xml``  | Galactic Plane Survey performed with threshold arrays  |
 +----------------------------+--------------------------------------------------------+
 | ``obs_egal_baseline.xml``  | Extragalactic Survey performed with baseline arrays    |
 +----------------------------+--------------------------------------------------------+
 | ``obs_egal_threshold.xml`` | Extragalactic Survey performed with threshold arrays   |
 +----------------------------+--------------------------------------------------------+
 | ``obs_all_baseline.xml``   | All surveys combined performed with baseline arrays    |
 +----------------------------+--------------------------------------------------------+
 | ``obs_all_threshold.xml``  | All surveys combined performed with threshold arrays   |
 +----------------------------+--------------------------------------------------------+


Models
^^^^^^

The ``models`` folder contains the definitions of all source and background
models that were used for simulating the data.
The file ``models_gc.xml`` is a so called
:ref:`Model Definition File <glossary_moddef>`
that collects the definition of all model components used for the Galactic
Centre Survey simulation.
The files ``models_gps.xml`` and ``models_egal.xml`` are the equivalent files
for the Galactic Plane Survey and the Extragalactic Survey.
The other files in the folder are FITS files containing sky maps that were
used as spatial templates for modelling extended or diffuse emission components.

.. warning::
   The FITS template files should always reside in the same folder as the
   :ref:`model definition XML files <glossary_moddef>`
   since the latter reference the former.



