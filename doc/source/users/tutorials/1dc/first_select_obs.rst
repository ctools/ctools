.. _1dc_first_select_obs:

Selecting the relevant observations
-----------------------------------

Let us in the following assume that you want to analyse the Galactic Plane
Survey data obtained with baseline arrays around the Galactic Centre.
For the purpose of this first analysis, you should create a subfolder under
the ``1dc`` folder (or any other convenient place) that will contain the
results for your first analysis.
You do this by typing:

.. code-block:: bash

   $ cd 1dc
   $ mkdir my_first_analysis
   $ cd my_first_analysis

The first step of the analysis consists in selecting the observations from the
:ref:`Observation Definition File <glossary_obsdef>`
``obs_gps_baseline.xml`` that have pointing directions close to the Galactic
Centre.
You do this with the :ref:`csobsselect` script by typing:

.. code-block:: bash

   $ csobsselect
   Input event list or observation definition XML file [obs.xml] ../data/obs_gps_baseline.xml
   Pointing selection region shape (CIRCLE|BOX) [CIRCLE]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] GAL
   Galactic longitude of selection centre (deg) (0-360) [184.56] 0.0
   Galactic latitude of selection centre (deg) (-90-90) [-5.79] 0.0
   Radius of selection circle (deg) (0-180) [5.0]
   Output observation definition XML file [outobs.xml] obs.xml

You selected ``CIRCLE`` as the shape of the pointing selection region,
specified the Galactic Centre as the centre of the selection region, and
selected all observations with pointing directions within 5 degrees from
the Galactic Centre from the ``obs_gps_baseline.xml`` file.
The selected observations are written into the file ``obs.xml`` that will be
created in the ``my_first_analysis`` folder.
The content of the file will look as follows:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="Galactic Plane Survey KSP" id="000020" instrument="CTA">
       <parameter name="EventList" file="../data/gps_baseline_000020.fits" />
       <parameter name="Calibration" database="prod2" response="South_50h" />
     </observation>
     <observation name="Galactic Plane Survey KSP" id="000021" instrument="CTA">
       <parameter name="EventList" file="../data/gps_baseline_000021.fits" />
       <parameter name="Calibration" database="prod2" response="South_50h" />
     </observation>
     <observation name="Galactic Plane Survey KSP" id="000022" instrument="CTA">
       <parameter name="EventList" file="../data/gps_baseline_000022.fits" />
       <parameter name="Calibration" database="prod2" response="South_50h" />
     </observation>
     ...
   </observation_list>

Each ``<observation>`` element corresponds to one observation that is identified
by a name attribute and a unique identifier attribute.
An ``<observation>`` element contains two parameters:
the ``EventList`` parameter that specifies the name of the corresponding event
file relative to the current working directory, and
the ``Calibration`` parameter that specifies the
:ref:`Instrument Response Function <glossary_irf>` that applies to the
event file.

