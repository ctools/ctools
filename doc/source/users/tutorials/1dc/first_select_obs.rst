.. _1dc_first_select_obs:

Selecting the relevant observations
-----------------------------------

Let's assume that you want to analyse the central region of our Galaxy using
the data obtained during the Galactic Centre Survey.

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
   Input event list or observation definition XML file [obs.xml] $CTADATA/obs/obs_gc_baseline.xml
   Pointing selection region shape (CIRCLE|BOX) [CIRCLE]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] GAL
   Galactic longitude of selection centre (deg) (0-360) [184.56] 0.0
   Galactic latitude of selection centre (deg) (-90-90) [-5.79] 0.0
   Radius of selection circle (deg) (0-180) [5.0] 0.1
   Output observation definition XML file [outobs.xml] obs.xml

You selected ``CIRCLE`` as the shape of the pointing selection region,
specified the Galactic Centre as the centre of the selection region, and
selected all observations with pointing directions within 0.1 degrees from
the Galactic Centre from the ``obs_gc_baseline.xml`` file.
The selected observations are written into the file ``obs.xml`` that will be
created in the ``my_first_analysis`` folder.
The content of the file will look as follows:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="GC" id="000414" instrument="CTA">
       <parameter name="EventList" file="/Users/jurgen/analysis/cta/dc/1dc/validation/milano-noedisp/1dc.pre/data/baseline/gc/gc_baseline_000414.fits.gz" />
       <parameter name="Calibration" database="prod3b" response="South_z20_50h" />
     </observation>
     <observation name="GC" id="000415" instrument="CTA">
       <parameter name="EventList" file="/Users/jurgen/analysis/cta/dc/1dc/validation/milano-noedisp/1dc.pre/data/baseline/gc/gc_baseline_000415.fits.gz" />
       <parameter name="Calibration" database="prod3b" response="South_z20_50h" />
     </observation>
     <observation name="GC" id="000416" instrument="CTA">
       <parameter name="EventList" file="/Users/jurgen/analysis/cta/dc/1dc/validation/milano-noedisp/1dc.pre/data/baseline/gc/gc_baseline_000416.fits.gz" />
       <parameter name="Calibration" database="prod3b" response="South_z20_50h" />
     </observation>
     ...
     <observation name="GC" id="000581" instrument="CTA">
       <parameter name="EventList" file="/Users/jurgen/analysis/cta/dc/1dc/validation/milano-noedisp/1dc.pre/data/baseline/gc/gc_baseline_000581.fits.gz" />
       <parameter name="Calibration" database="prod3b" response="South_z20_50h" />
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

