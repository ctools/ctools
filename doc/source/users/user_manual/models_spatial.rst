.. _um_models_spatial:

Spatial source model components
-------------------------------

The following sections present the spatial model components that are available 
in ctools for gamma-ray sources.

.. note::
   Except of the ``DiffuseMapCube`` model, all spatial models are normalised
   so that when integrated over the sphere the result is unity.

Point source
^^^^^^^^^^^^

  The ``PointSource`` model describes a point source

  .. code-block:: xml

    <source name="Crab" type="PointSource">
      <spatialModel type="PointSource">
        <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
        <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>

  where

  * ``RA`` is the Right Ascension (degrees)
  * ``DEC`` is the Declination (degrees)

  Instead of ``RA`` and ``DEC`` one may use ``GLON`` and ``GLAT`` to give the
  position in Galactic coordinates

  .. code-block:: xml

    <source name="Crab" type="PointSource">
      <spatialModel type="PointSource">
        <parameter name="GLON" scale="1.0" value="184.5575" min="-360" max="360" free="1"/>
        <parameter name="GLAT" scale="1.0" value="-5.7844"  min="-90"  max="90"  free="1"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>

  where

  * ``GLON`` is the Galactic longitude (degrees)
  * ``GLAT`` is the Galactic latitude (degrees)

  .. note::
     Galactic coordinates using the ``GLON`` and ``GLAT`` parameters can be
     specified for all spatial models that require a sky coordinate. Galactic
     coordinates are transformed internally into celestial coordinates, and
     in case that a model is written into an XML file, the corresponding
     parameters will be replaced by the celestial coordinates ``RA`` and
     ``DEC``.

  .. note::
     For compatibility with the Fermi/LAT ScienceTools the model type
     ``PointSource`` can be replaced by ``SkyDirFunction``.


Radial source
^^^^^^^^^^^^^

RadialDisk
~~~~~~~~~~

  The ``RadialDisk`` model describes a uniform intensity distribution within
  a given radius

  .. code-block:: xml

     <source name="Crab" type="ExtendedSource">
       <spatialModel type="RadialDisk">
         <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
         <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
         <parameter name="Radius" scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  where

  * ``RA`` is the Right Ascension of the disk centre (degrees)
  * ``DEC`` is the Declination of the disk centre (degrees)
  * ``Radius`` is the disk radius (degrees)

RadialRing
~~~~~~~~~~

  The ``RadialRing`` model specifies a uniform intensity distribution within
  a circular ring

  .. code-block:: xml

     <source name="Crab" type="ExtendedSource">
       <spatialModel type="RadialRing">
         <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
         <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
         <parameter name="Radius" scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>
         <parameter name="Width"  scale="1.0" value="0.15"    min="0.01" max="10"  free="1"/>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  where

  * ``RA`` is the Right Ascension of the ring centre (degrees)
  * ``DEC`` is the Declination of the ring centre (degrees)
  * ``Radius`` is the inner ring radius (degrees)
  * ``Width`` is the ring width radius (degrees)

  .. note::
     Specifying the inner ring radius and ring width guarantees that both
     parameters are well defined. The ring outer radius is given by
     ``Radius+Width``.


RadialGaussian
~~~~~~~~~~~~~~

  The ``RadialGaussian`` model describes a Gaussian intensity distribution

  .. code-block:: xml

     <source name="Crab" type="ExtendedSource">
       <spatialModel type="RadialGaussian">
         <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
         <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
         <parameter name="Sigma" scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  and implements

  .. math::
     M_{\rm spatial}(\theta) = \frac{1}{2 \pi \sigma^2} \exp
                    \left(-\frac{1}{2}\frac{\theta^2}{\sigma^2} \right),

  where

  * ``RA`` is the Right Ascension of the Gaussian centre (degrees)
  * ``DEC`` is the Declination of the Gaussian centre (degrees)
  * :math:`\sigma` = ``Sigma`` (degrees)

RadialShell
~~~~~~~~~~~

  The ``RadialShell`` model describes a spherical shell projected on the sky

  .. code-block:: xml

     <source name="Crab" type="ExtendedSource">
       <spatialModel type="RadialShell">
         <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
         <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
         <parameter name="Radius" scale="1.0" value="0.30"    min="0.01" max="10"  free="1"/>
         <parameter name="Width"  scale="1.0" value="0.10"    min="0.01" max="10"  free="1"/>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  and implements

  .. math::
     M_{\rm spatial}(\theta) =  n_0 \left \{
     \begin{array}{l l}
        \displaystyle
        \sqrt{ \theta_{\rm out}^2 - \theta^2 } - \sqrt{ \theta_{\rm in}^2 - \theta^2 }
        & \mbox{if $\theta \le \theta_{\rm in}$} \\
        \\
       \displaystyle
        \sqrt{ \theta_{\rm out}^2 - \theta^2 }
        & \mbox{if $\theta_{\rm in} < \theta \le \theta_{\rm out}$} \\
        \\
       \displaystyle
       0 & \mbox{if $\theta > \theta_{\rm out}$}
     \end{array}
     \right .

  where

  * ``RA`` is the Right Ascension of the shell centre (degrees)
  * ``DEC`` is the Declination of the shell centre (degrees)
  * :math:`\theta_{\rm out}` = ``Radius`` + ``Width`` (degrees)
  * :math:`\theta_{\rm in}` = ``Radius`` (degrees)


Elliptical source
^^^^^^^^^^^^^^^^^

EllipticalDisk
~~~~~~~~~~~~~~

  The ``EllipticalDisk`` model describes a uniform intensity distribution within
  an elliptical circumference:

  .. code-block:: xml

     <source name="Crab" type="ExtendedSource">
       <spatialModel type="EllipticalDisk">
         <parameter name="RA"          scale="1.0" value="83.6331" min="-360"  max="360" free="1"/>
         <parameter name="DEC"         scale="1.0" value="22.0145" min="-90"   max="90"  free="1"/>
         <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
         <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
         <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  where

  * ``RA`` is the Right Ascension (degrees)
  * ``DEC`` is the Declination (degrees)
  * ``PA`` is the position angle, counted counterclockwise from North (degrees)
  * ``MinorRadius`` is the minor radius of the ellipse (degrees)
  * ``MajorRadius`` is the major radius of the ellipse (degrees)

EllipticalGaussian
~~~~~~~~~~~~~~~~~~

  The ``EllipticalGaussian`` model describes a Gaussian intensity distribution

  .. code-block:: xml

    <source name="Crab" type="ExtendedSource">
      <spatialModel type="EllipticalGaussian">
        <parameter name="RA"          scale="1.0" value="83.6331" min="-360"  max="360" free="1"/>
        <parameter name="DEC"         scale="1.0" value="22.0145" min="-90"   max="90"  free="1"/>
        <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
        <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
        <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>

  and implements

  .. math::
     M_{\rm spatial}(\theta, \phi) = \exp \left( -\frac{\theta^2}{2 r_\mathrm{eff}^2} \right),

  with

  .. math::
     r_\mathrm{eff} = \frac{ab} {\sqrt{\left( a \sin (\phi - \phi_0) \right)^2 +
                      \sqrt{\left( b \cos (\phi - \phi_0) \right)^2}}}

  where

  * ``RA`` is the Right Ascension (degrees)
  * ``DEC`` is the Declination (degrees)
  * ``PA`` is the position angle, counted counterclockwise from North (degrees)
  * :math:`a` = ``MinorRadius`` (degrees)
  * :math:`b` = ``MajorRadius`` (degrees)
  * :math:`\phi_0` is the position angle of the ellipse, counted counterclockwise
    from North
  * :math:`\phi` is the azimuth angle with respect to North.


Diffuse source
^^^^^^^^^^^^^^

DiffuseIsotropic
~~~~~~~~~~~~~~~~

  The ``DiffuseIsotropic`` model describes an isotropic intensity distribution

  .. code-block:: xml

     <source name="Crab" type="DiffuseSource">
       <spatialModel type="DiffuseIsotropic">
         <parameter name="Value" scale="1" value="1" min="1"  max="1" free="0"/>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  where

  * ``Value`` is isotropic intensity

  .. note::
     For compatibility with the Fermi/LAT ScienceTools the model type
     ``DiffuseIsotropic`` can be replaced by ``ConstantValue``.

DiffuseMap
~~~~~~~~~~

  The ``DiffuseMap`` model describes an arbitrary intensity distribution in
  form of a sky map

  .. code-block:: xml

     <source name="Crab" type="DiffuseSource">
       <spatialModel type="DiffuseMap" file="map.fits">
         <parameter name="Normalization" scale="1" value="1" min="0.001" max="1000.0" free="0"/>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  where

  * ``Normalization`` is a normalization value

  and the ``file`` attribute specifies a sky map FITS file name. If a file name
  without path is specified it is assumed that the FITS file resides in the same
  directory as the model definition XML file.

  .. note::
     For compatibility with the Fermi/LAT ScienceTools the model type
     ``DiffuseMap`` can be replaced by ``SpatialMap`` and the parameter
     ``Normalization`` can be replaced by ``Prefactor``.

  .. note::
     By default, the diffuse map is normalised so that

     .. math::
        \int_{\Omega} M_{\rm spatial}(p|E) \, d\Omega = 1

     which means that the units of the spatial model component are
     :math:`[M_{\rm spatial}] = {\rm sr}^{-1}`.
     To avoid the normalisation you may add the ``normalize="0"`` attribute
     to the spatial model tag.

     .. code-block:: xml

        <source name="Crab" type="DiffuseSource">
          <spatialModel type="DiffuseMap" file="map.fits" normalize="0">
            <parameter name="Normalization" scale="1" value="1" min="0.001" max="1000.0" free="0"/>
          </spatialModel>
          <spectrum type="...">
            ...
          </spectrum>
        </source>

     In that case, generally

     .. math::
        \int_{\Omega} M_{\rm spatial}(p|E) \, d\Omega \neq 1

     and the spectral component cannot be directly interpreted as a physical
     source intensity.

  .. _um_models_spatial_src_diffuse_cube:

DiffuseMapCube
~~~~~~~~~~~~~~

  The ``DiffuseMapCube`` model describes an arbitrary energy-dependent intensity
  distribution in form of a map cube

  .. code-block:: xml

     <source name="Crab" type="DiffuseSource">
       <spatialModel type="DiffuseMapCube" file="map_cube.fits">
         <parameter name="Normalization" scale="1" value="1" min="0.001" max="1000.0" free="0"/>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  where

  * ``Normalization`` is a normalization value

  Note that the map cube is not normalised to unit, hence generally

  .. math::
     \int_{\Omega} M_{\rm spatial}(p|E) \, d\Omega \neq 1

  To compute the flux in a given energy band for a ``DiffuseMapCube`` model
  you have to integrated both the spatial and spectral components, i.e.

  .. math::
     \Phi = \int_{\Omega} \int_{E} M_{\rm spatial}(p|E) \times
            M_{\rm spectral}(E)\, d\Omega \, dE

  .. note::
     For compatibility with the Fermi/LAT ScienceTools the model type
     ``DiffuseMapCube`` can be replaced by ``MapCubeFunction`` and the parameter
     ``Normalization`` can be replaced by ``Value``.


Composite model
^^^^^^^^^^^^^^^

  The ``Composite`` model implements a composite model that is the sum of an
  arbitrary number of spatial models

  .. code-block:: xml

     <source name="Crab" type="CompositeSource">
       <spatialModel type="Composite">
         <spatialModel type="PointSource" component="PointSource">
           <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
           <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
         </spatialModel>
         <spatialModel type="RadialGaussian">
           <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
           <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
           <parameter name="Sigma" scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>
         </spatialModel>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  which implements

  .. math::
     M_{\rm spatial}(p|E) = \frac{1}{N} \sum_{i=0}^{N-1} M_{\rm spatial}^{(i)}(p|E)

  where :math:`M_{\rm spatial}^{(i)}(p|E)` is any spatial model component
  (including another composite model), and :math:`N` is the number of
  model components that are combined.
