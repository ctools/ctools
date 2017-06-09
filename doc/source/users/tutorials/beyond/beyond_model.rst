.. _sec_connecting_model:

Connecting observations to specific models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the previous examples there was always a single background model
component to describe the residual particle background in the
various dataset.
This implies that the spatial and spectral shape of the background
distribution is assumed to be identical for all observations.
This is fine in a simulation, but for a real life situation this
assumption will probably not hold.

Let's go back to the case of a joint analysis of two observations,
specified by the following observation definition XML file:

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <observation_list title="observation library">
    <observation name="Crab" id="00001" instrument="CTA">
      <parameter name="EventList" file="events1.fits"/>
    </observation>
    <observation name="Crab" id="00002" instrument="CTA">
      <parameter name="EventList" file="events2.fits"/>
    </observation>
  </observation_list>

Each observation has a **unique identifier string** and you can use this 
unique string to connect an observation to a specific model component.
This is illustrated in the model XML file below that now has two background
components.
Both components have an ``id`` attribute that ties the model to a specific
observation.
Using ``id="00001"``, the first background model is tied to the 
``events1.fits`` file, using ``id="00002"``, the second background model is
tied to the ``events2.fits`` file.

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <source_library title="source library">
    <source name="Crab" type="PointSource">
      <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"       scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
      </spectrum>
      <spatialModel type="PointSource">
        <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="0"/>
        <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="0"/>
      </spatialModel>
    </source>
    <source name="Background_00001" type="CTAIrfBackground" instrument="CTA" id="00001">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor"   scale="1.0"  value="1.0"  min="1e-3" max="1e+3"   free="1"/>
        <parameter name="Index"       scale="1.0"  value="0.0"  min="-5.0" max="+5.0"   free="1"/>
        <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
      </spectrum>
    </source>
    <source name="Background_00002" type="CTAIrfBackground" instrument="CTA" id="00002">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor"   scale="1.0"  value="10.0"  min="1e-3" max="1e+3"   free="1"/>
        <parameter name="Index"       scale="1.0"  value="0.0"  min="-5.0" max="+5.0"   free="1"/>
        <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
      </spectrum>
    </source>
  </source_library>

.. note::

   **Model components in a XML file need to have a unique name.**
   For this reason the components are named here ``Background_00001``
   and ``Background_00002``.

Note that both background model components distinguish in their event rate.
``Background_00002`` has a ten times larger background rate than
``Background_00001``.
To illustrate this difference, two observations offset by +/- 2 degrees from
the Crab nebulae have been simulated using :ref:`ctobssim`.
A joint fit of the observations is then performed using :ref:`ctlike` and 
a residual map is created using :ref:`csresmap`.

The result is shown below.
Obviously, the upper observation has larger residuals owing to the ten 
times larger background rate that was used in the simulation of the
associated data.

.. figure:: resmap-bkg.png
   :height: 400px
   :align: center

   *Residual map of two jointly analysed observations offset by +/- 2 degrees*

The scheme is even more versatile in that it allows to connect a given 
model component to several specific observations.
This is illustrated in the model definition file below where component
``Background_00001`` is now applicable for observations ``00001``, 
``00003`` and ``00004``.

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <source_library title="source library">
    <source name="Crab" type="PointSource">
      <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"       scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
      </spectrum>
      <spatialModel type="PointSource">
        <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="0"/>
        <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="0"/>
      </spatialModel>
    </source>
    <source name="Background_00001" type="CTAIrfBackground" instrument="CTA" id="00001,00003,00004">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor"   scale="1.0"  value="1.0"  min="1e-3" max="1e+3"   free="1"/>
        <parameter name="Index"       scale="1.0"  value="0.0"  min="-5.0" max="+5.0"   free="1"/>
        <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
      </spectrum>
    </source>
    <source name="Background_00002" type="CTAIrfBackground" instrument="CTA" id="00002">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor"   scale="1.0"  value="10.0"  min="1e-3" max="1e+3"   free="1"/>
        <parameter name="Index"       scale="1.0"  value="0.0"  min="-5.0" max="+5.0"   free="1"/>
        <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
      </spectrum>
    </source>
  </source_library>



   

