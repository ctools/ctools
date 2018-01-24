.. _howto_extent:

How to determine the extension of a source?
-------------------------------------------

  .. admonition:: What you will learn

     You will learn how you **determine the spatial extent of a source** by
     fitting an extended source model to the data.

To determine the spatial extension of a source, an extended source model
has to be used in the
:ref:`model definition file <glossary_moddef>`
for the source of interest.
The example below is based on the source model that was determined
during the
:ref:`iterative model improvement <1dc_first_improving>`
of the
:ref:`first CTA Data Challenge <glossary_1dc>`
tutorial where two point sources were detected.
As in the tutorial, ``Src001`` is modelled with an exponentially cut-off power
law spectrum, ``Src002`` with a simple power law.
In addition, a diffuse emission component was added to the model.

As noted during the
:ref:`iterative improvement of the source model <1dc_first_improving>`
there is a bright extended source south-west of the Galactic Centre, and
we add now an additional ``Src003`` component with a ``RadialDisk`` model to
describe the spatial morphology of the source.
The
:ref:`model definition <glossary_moddef>`
for ``Src003`` is shown below:

.. code-block:: xml

   <source name="Src003" type="ExtendedSource" tscalc="1">
     <spectrum type="PowerLaw">
       <parameter name="Prefactor"   scale="5.7e-18" value="1.0" min="0"    max="1000.0" free="1"/>
       <parameter name="Index"       scale="-2.48"   value="1.0" min="-4.0" max="4.0"    free="1"/>
       <parameter name="PivotEnergy" scale="300000"  value="1.0" free="0" />
     </spectrum>
     <spatialModel type="RadialDisk">
       <parameter name="RA"     scale="1.0" value="266.3070" min="-360"  max="360" free="1"/>
       <parameter name="DEC"    scale="1.0" value="-30.1876" min="-90"   max="90"  free="1"/>
       <parameter name="Radius" scale="1.0" value="0.3"      min="0.001" max="10"  free="1"/>
     </spatialModel>
   </source>

Running :ref:`ctlike` using the
:ref:`model definition file <glossary_moddef>`
will then fit the disk radius of ``Src003``. The results can be seen in the
log file and will also be written into the output
:ref:`model definition file <glossary_moddef>`.
For illustration, an excerpt of the :ref:`ctlike` log file is shown below.
The disk with of ``Src003`` has been fitted to 0.335 +/- 0.003 degrees.
The source is detected with a TS value of 5393.8.

.. code-block:: none

   2017-08-25T19:37:47: === GModelSky ===
   2017-08-25T19:37:47:  Name ......................: Src003
   2017-08-25T19:37:47:  Instruments ...............: all
   2017-08-25T19:37:47:  Test Statistic ............: 5393.78645148844
   2017-08-25T19:37:47:  Instrument scale factors ..: unity
   2017-08-25T19:37:47:  Observation identifiers ...: all
   2017-08-25T19:37:47:  Model type ................: ExtendedSource
   2017-08-25T19:37:47:  Model components ..........: "RadialDisk" * "PowerLaw" * "Constant"
   2017-08-25T19:37:47:  Number of parameters ......: 7
   2017-08-25T19:37:47:  Number of spatial par's ...: 3
   2017-08-25T19:37:47:   RA .......................: 266.306975228126 +/- 0.00363322258453142 [-360,360] deg (free,scale=1)
   2017-08-25T19:37:47:   DEC ......................: -30.1946569442618 +/- 0.0031545381615373 [-90,90] deg (free,scale=1)
   2017-08-25T19:37:47:   Radius ...................: 0.335173109280516 +/- 0.00257226134023012 [0.001,10] deg (free,scale=1)
   2017-08-25T19:37:47:  Number of spectral par's ..: 3
   2017-08-25T19:37:47:   Prefactor ................: 6.23930379041412e-17 +/- 1.06204679984262e-18 [0,5.7e-15] ph/cm2/s/MeV (free,scale=5.7e-18,gradient)
   2017-08-25T19:37:47:   Index ....................: -2.73810250489576 +/- 0.0151068737641297 [9.92,-9.92]  (free,scale=-2.48,gradient)
   2017-08-25T19:37:47:   PivotEnergy ..............: 300000 MeV (fixed,scale=300000,gradient)
   2017-08-25T19:37:47:  Number of temporal par's ..: 1
   2017-08-25T19:37:47:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

.. warning::
   The parameter errors returned by :ref:`ctlike` are purely statistical. In
   addition to the statistical errors there are systematic uncertainties, such
   as for example the limits on the knowledge of the point spread function.
   These systematic uncertainties are **not** determined by ctools.

The figure below shows the residual map, generated using the ``SIGNIFICANCE``
method, after subtracting the three fitted
sources and the diffuse emission model from the data. There is a ring-like
residual at the position of ``Src003`` which suggests that an axisymmetric
disk is not an accurate description of the data.

.. figure:: howto_extent_disk.png
   :width: 400px
   :align: center

   *Residual sky map for a radial disk spatial shape for Src003*

There are other spatial models in ctools, and we try in a second step an
axisymmetric Gaussian spatial shape. The corresponding
:ref:`model definition file <glossary_moddef>`
looks as follows:

.. code-block:: xml

   <source name="Src003" type="ExtendedSource" tscalc="1">
     <spectrum type="PowerLaw">
       <parameter name="Prefactor"   scale="5.7e-18" value="1.0" min="0"    max="1000.0" free="1"/>
       <parameter name="Index"       scale="-2.48"   value="1.0" min="-4.0" max="4.0"    free="1"/>
       <parameter name="PivotEnergy" scale="300000"  value="1.0" free="0" />
     </spectrum>
     <spatialModel type="RadialGaussian">
       <parameter name="RA"    scale="1.0" value="266.3070" min="-360"  max="360" free="1"/>
       <parameter name="DEC"   scale="1.0" value="-30.1876" min="-90"   max="90"  free="1"/>
       <parameter name="Sigma" scale="1.0" value="0.2"      min="0.001" max="10"  free="1"/>
     </spatialModel>
   </source>

Running :ref:`ctlike` again with that model results in a Gaussian sigma of
0.198 +/- 0.003 degrees for ``Src003``. The source is detected with a TS value
of 5665.7 which is considerably larger than the TS value of 5393.8 that is found
above for the radial disk model.
Below an excerpt of the :ref:`ctlike` log file:

.. code-block:: none

   2017-08-25T20:32:55: === GModelSky ===
   2017-08-25T20:32:55:  Name ......................: Src003
   2017-08-25T20:32:55:  Instruments ...............: all
   2017-08-25T20:32:55:  Test Statistic ............: 5665.73173788155
   2017-08-25T20:32:55:  Instrument scale factors ..: unity
   2017-08-25T20:32:55:  Observation identifiers ...: all
   2017-08-25T20:32:55:  Model type ................: ExtendedSource
   2017-08-25T20:32:55:  Model components ..........: "RadialGaussian" * "PowerLaw" * "Constant"
   2017-08-25T20:32:55:  Number of parameters ......: 7
   2017-08-25T20:32:55:  Number of spatial par's ...: 3
   2017-08-25T20:32:55:   RA .......................: 266.30022317901 +/- 0.00475655134090925 [-360,360] deg (free,scale=1)
   2017-08-25T20:32:55:   DEC ......................: -30.1993627725406 +/- 0.0041251924859714 [-90,90] deg (free,scale=1)
   2017-08-25T20:32:55:   Sigma ....................: 0.198343910095786 +/- 0.00298438230006778 [0.001,10] deg (free,scale=1)
   2017-08-25T20:32:55:  Number of spectral par's ..: 3
   2017-08-25T20:32:55:   Prefactor ................: 7.15427332759485e-17 +/- 1.40842792023593e-18 [0,5.7e-15] ph/cm2/s/MeV (free,scale=5.7e-18,gradient)
   2017-08-25T20:32:55:   Index ....................: -2.71899753120954 +/- 0.0146676775081257 [9.92,-9.92]  (free,scale=-2.48,gradient)
   2017-08-25T20:32:55:   PivotEnergy ..............: 300000 MeV (fixed,scale=300000,gradient)
   2017-08-25T20:32:55:  Number of temporal par's ..: 1
   2017-08-25T20:32:55:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

.. note::
   While the TS values can formally not be convert into a statistical
   significance between different spatial model hypotheses, a TS improvement
   of 271.9 indicates a considerably better fit of the axisymmetric Gaussian
   model with respect to the radial disk model to the data.

The figure below shows the residual map for the fit of ``Src003`` with an
axisymmetric Gaussian model.
The map now looks pretty flat around ``Src003``, suggesting that an axisymmetric
Gaussian model is an appropriate description for the morphology of the gamma-ray
source.

.. figure:: howto_extent_gauss.png
   :width: 400px
   :align: center

   *Residual sky map for an axisymmetric Gaussian spatial shape for Src003*

.. tip::
   The region overlays for the residual sky maps were generated using the
   :ref:`csmodelinfo` script.

.. warning::
   The fitting of extended spatial models takes more computing time
   than the fitting of point sources. The computing time is related to the
   spatial extent of the source and to the spatial shape, with a Gaussian
   disk model taking considerably more computing time than a radial disk
   model due to the tails of the Gaussian function. It is therefore **recommended
   to use by default radial disk models for the extension fitting**, and only
   switch to a Gaussian disk models when really needed, or for the determination
   of final values for a publication.

   For reference, here the computing times on Mac OS X for the example using
   different spatial morphology hypotheses for ``Src003``:

   * Point source: 12 min
   * Radial disk source: 15 min
   * Axisymmetric Gaussian source: 40 min
