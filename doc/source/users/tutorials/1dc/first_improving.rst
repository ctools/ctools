.. _1dc_first_improving:

Iteratively improving the model
-------------------------------

  .. admonition:: What you will learn

     You will learn how to **iteratively improve** your source model by
     inspecting the fit residuals and by adjusting the model components as needed.

After doing the model fit you should investigate the residuals to verify that
the model components properly describe the observed event distribution.
You do this with the :ref:`csresmap` script by providing the output
:ref:`model definition file <glossary_moddef>`
``results.xml`` produced by :ref:`ctlike` as input:

.. code-block:: bash

   $ csresmap
   Input event list, counts cube, or observation definition XML file [events.fits] cntcube.fits
   Input model cube file (generated with ctmodel) [expcube.fits] NONE
   Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits
   Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits
   Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits
   Input model definition XML file [$CTOOLS/share/models/crab.xml] stacked_results_cutoff.xml
   Residual map computation algorithm (SUB|SUBDIV|SUBDIVSQRT|SIGNIFICANCE) [SUBDIV] SUB
   Output residual map file [resmap.fits] resmap.fits

This produces the file ``resmap.fits`` that contains a residual map that
you can display for example with `ds9 <http://ds9.si.edu>`_.
The figure below shows the map, with the fitted source positions overlayed
as white circles.

.. figure:: first_skymap_residual.png
   :width: 400px
   :align: center

   *Residual sky map after subtraction of the fitted model*

The subtraction of the fitted sources clearly leaves some emission holes at
their locations, which means that the fitted sources picked up some of the
underlying Galactic diffuse gamma-ray emission. To prevent this cross-talk,
a diffuse emission model should be added to the
:ref:`model definition file <glossary_moddef>`,
as illustrated below:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Src001" type="PointSource">
       <spectrum type="ExponentialCutoffPowerLaw">
         <parameter name="Prefactor"    scale="1e-18" value="5.7"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"        scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="CutoffEnergy" scale="1e6"   value="10.0" min="0.01"  max="1000.0" free="1"/>
         <parameter name="PivotEnergy"  scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA" value="266.4120906928" error="0" scale="1" free="1" />
         <parameter name="DEC" value="-29.0219729468991" error="0" scale="1" free="1" />
       </spatialModel>
     </source>
     <source name="Src002" type="PointSource">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" value="1" error="0" scale="5.7e-18" min="0" free="1" />
         <parameter name="Index" value="1" error="-0" scale="-2.48" min="-4.03225806451613" max="4.03225806451613" free="1" />
         <parameter name="PivotEnergy" value="1" scale="300000" free="0" />
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA" value="266.793148201606" error="0" scale="1" free="1" />
         <parameter name="DEC" value="-28.1253038992376" error="0" scale="1" free="1" />
       </spatialModel>
     </source>
     <source name="IEM" type="DiffuseSource">
       <spectrum type="ConstantValue">
         <parameter name="Value" value="1" scale="1" min="1e-05" max="100000" free="1" />
       </spectrum>
       <spatialModel type="MapCubeFunction" file="$CTADATA/models/cube_iem.fits.gz">
         <parameter name="Normalization" value="1" scale="1" min="0.001" max="1000" free="0" />
       </spatialModel>
     </source>
     <source name="Background" type="CTAIrfBackground">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   value="1" scale="1" min="0.1" max="10" free="1" />
         <parameter name="Index"       value="0" scale="1" min="-10" max="10" free="1" />
         <parameter name="PivotEnergy" value="1" scale="1000000" free="0" />
       </spectrum>
     </source>
   </source_library>

Repeating the fit with this model and producing a corresponding residual map
produces the map shown below. Now, the residuals near the two point sources
are flatter and the diffuse emission has disappeared. However, there
is still a clear emission hole near ``Src001`` which suggests that this source
may be extended.

.. figure:: first_skymap_residual_iem.png
   :width: 400px
   :align: center

   *Residual sky map after subtraction of the fitted model including a diffuse emission component*

In our attempt to iteratively refine the source model, we therefore use in our next
iteration a radial disk source instead of a point source for ``Src001`` by
specifying the following
:ref:`model definition file <glossary_moddef>`
on input to the :ref:`ctlike` tool:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Src001" type="PointSource">
       <spectrum type="ExponentialCutoffPowerLaw">
         <parameter name="Prefactor"    scale="1e-18" value="5.7"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"        scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="CutoffEnergy" scale="1e6"   value="10.0" min="0.01"  max="1000.0" free="1"/>
         <parameter name="PivotEnergy"  scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
       </spectrum>
       <spatialModel type="RadialDisk">
         <parameter name="RA"     value="266.4121" scale="1.0" min="-360"  max="360" free="1"/>
         <parameter name="DEC"    value="-29.0220" scale="1.0" min="-90"   max="90"  free="1"/>
         <parameter name="Radius" value="0.01"     scale="1.0" min="0.001" max="10"  free="1"/>
       </spatialModel>
     </source>
     <source name="Src002" type="PointSource">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" value="1" error="0" scale="5.7e-18" min="0" free="1" />
         <parameter name="Index" value="1" error="-0" scale="-2.48" min="-4.03225806451613" max="4.03225806451613" free="1" />
         <parameter name="PivotEnergy" value="1" scale="300000" free="0" />
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA" value="266.793148201606" error="0" scale="1" free="1" />
         <parameter name="DEC" value="-28.1253038992376" error="0" scale="1" free="1" />
       </spatialModel>
     </source>
     <source name="IEM" type="DiffuseSource">
       <spectrum type="ConstantValue">
         <parameter name="Value" value="1" scale="1" min="1e-05" max="100000" free="1" />
       </spectrum>
       <spatialModel type="MapCubeFunction" file="$CTADATA/models/cube_iem.fits.gz">
         <parameter name="Normalization" value="1" scale="1" min="0.001" max="1000" free="0" />
       </spatialModel>
     </source>
     <source name="Background" type="CTAIrfBackground">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   value="1" scale="1" min="0.1" max="10" free="1" />
         <parameter name="Index"       value="0" scale="1" min="-10" max="10" free="1" />
         <parameter name="PivotEnergy" value="1" scale="1000000" free="0" />
       </spectrum>
     </source>
   </source_library>

As shown in the residual map below, the emission hole at the position of
``Src001`` has now disappeared. The disk model has in fact been fitted with
a radius of 5.3 +/- 0.1 arcmin, indicating that ``Src001`` is significantly
extended. The remaining residuals suggest that also ``Src002`` may be
slightly extended, but it is left as an exercise to the User to test this
hypothesis.

.. figure:: first_skymap_residual_iem_disk.png
   :width: 400px
   :align: center

   *Residual sky map after subtraction of the fitted model including a diffuse emission component and a disk model for Src001*
