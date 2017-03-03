.. _1dc_howto_extent:

How to determine the extension of a source?
-------------------------------------------

To determine the spatial extension of a source, an extended source model
has to be used in the
:ref:`model definition file <glossary_moddef>`
of the source of interest.
In the example below, a ``RadialGaussian`` model is used for the source
``Src001`` to model the spatial morphology of the source by an axisymmetric
Gaussian profile with the Gaussian width ``sigma`` as free parameter.
Note that at the same time the source type has been set to
``ExtendedSource``.

.. code-block:: xml

   <source name="Src001" type="ExtendedSource">
     <spectrum type="PowerLaw">
       <parameter name="Prefactor"   value="1" scale="5.7e-18" min="0" free="1" />
       <parameter name="Index"       value="1" scale="-2.48"   min="-5" max="5" free="1" />
       <parameter name="PivotEnergy" value="1" scale="300000" free="0" />
     </spectrum>
     <spatialModel type="RadialGaussian">
       <parameter name="RA"    scale="1.0" value="266.4240" min="-360"  max="360" free="1"/>
       <parameter name="DEC"   scale="1.0" value="-29.0049" min="-90"   max="90"  free="1"/>
       <parameter name="Sigma" scale="1.0" value="0.05"     min="0.001" max="10"  free="1"/>
     </spatialModel>
   </source>

Running :ref:`ctlike` using the modified
:ref:`model definition file <glossary_moddef>`
will then fit the Gaussian width of ``Src001``. The results can be seen in the
log file and will also be written into the output
:ref:`model definition file <glossary_moddef>`.

For illustration, an excerpt of the :ref:`ctlike` log file is shown below.
The Gaussian width of ``Src001`` has been fitted to 0.0218 +/- 0.0006 degrees.

.. code-block:: xml

   2017-03-04T17:44:06: === GModelSky ===
   2017-03-04T17:44:06:  Name ......................: Src001
   2017-03-04T17:44:06:  Instruments ...............: all
   2017-03-04T17:44:06:  Instrument scale factors ..: unity
   2017-03-04T17:44:06:  Observation identifiers ...: all
   2017-03-04T17:44:06:  Model type ................: ExtendedSource
   2017-03-04T17:44:06:  Model components ..........: "RadialGaussian" * "PowerLaw" * "Constant"
   2017-03-04T17:44:06:  Number of parameters ......: 7
   2017-03-04T17:44:06:  Number of spatial par's ...: 3
   2017-03-04T17:44:06:   RA .......................: 266.416716505268 +/- 0.000753191747203773 [-360,360] deg (free,scale=1)
   2017-03-04T17:44:06:   DEC ......................: -29.0044104159472 +/- 0.000655420291471696 [-90,90] deg (free,scale=1)
   2017-03-04T17:44:06:   Sigma ....................: 0.02183871432113 +/- 0.000636651479386333 [0.001,10] deg (free,scale=1)
   2017-03-04T17:44:06:  Number of spectral par's ..: 3
   2017-03-04T17:44:06:   Prefactor ................: 4.85894838349519e-17 +/- 7.40461052122477e-19 [0,infty[ ph/cm2/s/MeV (free,scale=5.7e-18,gradient)
   2017-03-04T17:44:06:   Index ....................: -2.29168205335762 +/- 0.00909936860902909 [10,-10]  (free,scale=-2.48,gradient)
   2017-03-04T17:44:06:   PivotEnergy ..............: 300000 MeV (fixed,scale=300000,gradient)
   2017-03-04T17:44:06:  Number of temporal par's ..: 1
   2017-03-04T17:44:06:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

.. warning::
   The errors returned by :ref:`ctlike` are purely statistical. In addition
   there are systematic uncertainties, such as for example the limits in the
   knowledge of the point spread function. It remains to be seen whether CTA
   will be able to reliably determine source extents as small as 0.02 degrees.
