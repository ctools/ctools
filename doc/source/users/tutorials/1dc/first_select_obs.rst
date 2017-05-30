.. _1dc_first_select_obs:

Selecting the relevant observations
-----------------------------------

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
:ref:`Observation Definition File <glossary_obsdef>`
``obs_gc_baseline.xml`` that have pointing directions close to the Galactic
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
The file will contain the definition of 57 observations and will look as
follows:

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
:ref:`Instrument Response Function <glossary_irf>` that applies to the
event file.

