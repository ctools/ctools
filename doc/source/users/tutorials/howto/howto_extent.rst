.. _1dc_howto_extent:

How to determine the extension of a source?
-------------------------------------------

  .. admonition:: What you will learn

     You will learn how you **determine the spatial extent of a source** by
     fitting an extended source model to the data.

To determine the spatial extension of a source, an extended source model
has to be used in the
:ref:`model definition file <glossary_moddef>`
of the source of interest.
In the example below, a ``RadialGaussian`` model is used for the source
``Src001`` to model the spatial morphology of the source by an axisymmetric
Gaussian profile with the Gaussian width ``sigma`` as free parameter.
Note that at the same time the source type was set to
``ExtendedSource``.

.. code-block:: xml

   <source name="Src001" type="PointSource">
     <spectrum type="ExponentialCutoffPowerLaw">
       <parameter name="Prefactor"    scale="1e-18" value="5.7"  min="1e-07" max="1000.0" free="1"/>
       <parameter name="Index"        scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
       <parameter name="CutoffEnergy" scale="1e6"   value="10.0" min="0.01"  max="1000.0" free="1"/>
       <parameter name="PivotEnergy"  scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
     </spectrum>
     <spatialModel type="RadialGaussian">
       <parameter name="RA"    scale="1.0" value="266.4121" min="-360"  max="360" free="1"/>
       <parameter name="DEC"   scale="1.0" value="-29.0220" min="-90"   max="90"  free="1"/>
       <parameter name="Sigma" scale="1.0" value="0.01"     min="0.001" max="10"  free="1"/>
     </spatialModel>
   </source>

Running :ref:`ctlike` using the modified
:ref:`model definition file <glossary_moddef>`
will then fit the Gaussian width of ``Src001``. The results can be seen in the
log file and will also be written into the output
:ref:`model definition file <glossary_moddef>`.

For illustration, an excerpt of the :ref:`ctlike` log file is shown below.
The Gaussian width of ``Src001`` has been fitted to 0.048 +/- 0.001 degrees.

.. code-block:: xml

   2017-06-08T15:07:46: === GModelSky ===
   2017-06-08T15:07:46:  Name ......................: Src001
   2017-06-08T15:07:46:  Instruments ...............: all
   2017-06-08T15:07:46:  Instrument scale factors ..: unity
   2017-06-08T15:07:46:  Observation identifiers ...: all
   2017-06-08T15:07:46:  Model type ................: ExtendedSource
   2017-06-08T15:07:46:  Model components ..........: "RadialGaussian" * "ExponentialCutoffPowerLaw" * "Constant"
   2017-06-08T15:07:46:  Number of parameters ......: 8
   2017-06-08T15:07:46:  Number of spatial par's ...: 3
   2017-06-08T15:07:46:   RA .......................: 266.417513130306 +/- 0.00132553002208653 [-360,360] deg (free,scale=1)
   2017-06-08T15:07:46:   DEC ......................: -29.0076435535415 +/- 0.00117490755596776 [-90,90] deg (free,scale=1)
   2017-06-08T15:07:46:   Sigma ....................: 0.0484844490479065 +/- 0.000997962295674141 [0.001,10] deg (free,scale=1)
   2017-06-08T15:07:46:  Number of spectral par's ..: 4
   2017-06-08T15:07:46:   Prefactor ................: 3.26956513671326e-17 +/- 6.26524553994195e-19 [1e-25,1e-15] ph/cm2/s/MeV (free,scale=1e-18,gradient)
   2017-06-08T15:07:46:   Index ....................: -2.32005423294123 +/- 0.0173504097487491 [-0,-5]  (free,scale=-1,gradient)
   2017-06-08T15:07:46:   CutoffEnergy .............: 162314916.632706 +/- 182782066.654979 [10000,1000000000] MeV (free,scale=1000000,gradient)
   2017-06-08T15:07:46:   PivotEnergy ..............: 300000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2017-06-08T15:07:46:  Number of temporal par's ..: 1
   2017-06-08T15:07:46:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

.. warning::
   The parameter errors returned by :ref:`ctlike` are purely statistical. In
   addition to the statistical errors there are systematic uncertainties, such
   as for example the limits on the knowledge of the point spread function.


