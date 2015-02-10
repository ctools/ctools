.. _models:

Modelling CTA data
------------------

.. _sec_model_concept:

Overall concept
~~~~~~~~~~~~~~~

A central element of any data analysis is the modeling of the observed 
event distribution.
Generally, the measured events can be separated into two distinct classes:
events being attributed to gamma rays of celestial origin (the "source" 
events) or events being attributed to any other kind of trigger (the 
"background" events).
In the first case, the :ref:`Instrument Response Functions (IRFs) <response>`
describe how an incident gamma-ray converts into a measured event.
In the second case there is often no general prescription, and the 
distribution of background events is commonly modeled directly in the data 
space of the observable quantities.
The main difficulty of gamma-ray astronomy is that not all events can be
tagged a priori as being source or background, hence their separation can 
only be done on a statistical basis, based on the different morphologies, 
spectral characteristics or eventually temporal signatures of both event
categories.

For this purpose, ctools use a general model that describes the spatial, 
spectral and temporal properties of the source and background components.
The model is composed of an arbitrary number of model components that
add up linearly to provide a prediction of the expected number of events.
Model components are generally parametric, and model parameters can be 
adjusted through a maximum likelihood procedure to find the set of 
parameters that represent best the measured data.
Model components representing celestial sources are convolved with the 
:ref:`IRFs <response>` to predict the expected number of source events in 
the data.
Background model components will be directly expressed as expected number 
of background events without any :ref:`IRF <response>` convolution.


.. _sec_model_implementation:

Implementation
~~~~~~~~~~~~~~

The general model is describe in ctools using a model definition XML file. 
Below is a simple example of such a file comprising one source and one 
background model.
Each model is factorised in a spectral (tag ``<spectrum>``) and a 
spatial component (tags ``<spatialModel>`` and ``<radialModel>``):

.. math::
  M(x,y,E) = M_{\rm spectral}(E) \times M_{\rm spatial}(x,y)

In this specific example, the source component ``Crab`` describes 
a point source at the location of the Crab nebula with a power law spectral
shape.
The background component ``Background`` is modelled as a radial Gaussian 
function in offset angle squared (with the offset angle being defined as 
the angle between pointing and measured event direction) and a spectral
function that is tabulated in an ASCII file.

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <source_library title="source library">
    <source name="Crab" type="PointSource">
      <spectrum type="PowerLaw">
         <parameter name="Prefactor" scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"     scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="Scale"     scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
      </spectrum>
      <spatialModel type="SkyDirFunction">
        <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
        <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
      </spatialModel>
    </source>
    <source name="Background" type="RadialAcceptance" instrument="CTA">
      <spectrum type="FileFunction" file="$CTOOLS/share/models/bkg_dummy.txt">
        <parameter name="Normalization" scale="1.0" value="1.0" min="0.0" max="1000.0" free="1"/>
      </spectrum>
      <radialModel type="Gaussian">
         <parameter name="Sigma" scale="1.0" value="3.0" min="0.01" max="10.0" free="1"/>
      </radialModel>
    </source>
  </source_library>

Model parameters are specified by a ``<parameter>`` tag with a certain 
number of attributes.
The ``name`` attribute gives the (case-sensitive) parameter name that 
needs to be unique within a given model component.
The ``scale`` attribute gives a scaling factor that will be multiplied by 
the ``value`` attribute to provide the real (physical) parameter value.
The ``min`` and ``max`` attributes specify boundaries for the ``value``
term of the parameter.
And the ``free`` attribute specifies whether a parameter should be fitted 
(``free="1"``) or kept fixed (``free="0"``) in a maximum likelihood 
analysis.
After a maximum likelihood fit using :ref:`ctlike`, an
``error`` attribute giving the statisical uncertainty of the ``value``
term will be appended to each ``<parameter>`` tag.

.. note::

   XML files are ASCII files and can be edited by hand using any text 
   editor.
   The indentation of the XML elements is not mandatory.

.. note::

   The splitting of parameter values into a ``value`` and ``scale`` factor 
   is mainly for numerical purposes.
   Parameter fitting algorithms can be ill-conditioned if several 
   parameters of vastly different orders of magnitudes need to be 
   optimised simultaneously.
   Splitting a value into two components allows to "prescale" the 
   parameters so that the effective parameters to be optimised (the ``value`` terms) 
   are all of about unity.

.. warning::

   The syntax of the model definition XML file has been adopted from the 
   syntax used by the Fermi/LAT ScienceTools.
   This should allow for exchange of model definition files between ctools 
   and the Fermi/LAT ScienceTools.
   Note, however, that ctools implements a number of extensions with 
   respect to the Fermi/LAT ScienceTools syntax, hence eventually minor 
   manual adjustments may be needed to use ctools XML files for 
   Fermi/LAT analyses.


.. _sec_spatial_models:

Spatial model components
~~~~~~~~~~~~~~~~~~~~~~~~

The following sections present the spatial model components that are available 
in ctools.

Point source
^^^^^^^^^^^^

  .. code-block:: xml

    <source name="Crab" type="PointSource">
      <spatialModel type="SkyDirFunction">
        <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
        <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>


Radial source
^^^^^^^^^^^^^

  .. code-block:: xml

    <source name="Crab" type="ExtendedSource">
      <spatialModel type="DiskFunction">
        <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
        <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
        <parameter name="Radius" scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>

  .. code-block:: xml

    <source name="Crab" type="ExtendedSource">
      <spatialModel type="GaussFunction">
        <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
        <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
        <parameter name="Sigma" scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>

  .. code-block:: xml

    <source name="Crab" type="ExtendedSource">
      <spatialModel type="ShellFunction">
        <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
        <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
        <parameter name="Radius" scale="1.0" value="0.30"    min="0.01" max="10"  free="1"/>
        <parameter name="Width"  scale="1.0" value="0.10"    min="0.01" max="10"  free="1"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>


Elliptical source
^^^^^^^^^^^^^^^^^

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


Diffuse source
^^^^^^^^^^^^^^

  .. code-block:: xml

    <source name="Crab" type="DiffuseSource">
      <spatialModel type="ConstantValue">
         <parameter name="Value" scale="1" value="1" min="1"  max="1" free="0"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>

  .. code-block:: xml

    <source name="Crab" type="DiffuseSource">
      <spatialModel type="SpatialMap" file="map.fits">
         <parameter name="Prefactor" scale="1" value="1" min="0.001" max="1000.0" free="0"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>

  .. code-block:: xml

    <source name="Crab" type="DiffuseSource">
      <spatialModel type="MapCubeFunction" file="map_cube.fits">
        <parameter name="Normalization" scale="1" value="1" min="0.001" max="1000.0" free="0"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>


CTA radial background
^^^^^^^^^^^^^^^^^^^^^

  .. code-block:: xml

    <source name="Background" type="RadialAcceptance" instrument="CTA">
      <radialModel type="Gaussian">
        <parameter name="Sigma" scale="1.0" value="3.0" min="0.01" max="10.0" free="1"/>
      </radialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>

  .. code-block:: xml

    <source name="Background" type="RadialAcceptance" instrument="CTA">
      <radialModel type="Profile">
        <parameter name="Width" scale="1.0" value="1.5" min="0.1" max="1000.0" free="1"/>
        <parameter name="Core"  scale="1.0" value="3.0" min="0.1" max="1000.0" free="1"/>
        <parameter name="Tail"  scale="1.0" value="5.0" min="0.1" max="1000.0" free="1"/>
      </radialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>

  .. code-block:: xml

    <source name="Background" type="RadialAcceptance" instrument="CTA">
      <radialModel type="Polynom">
        <parameter name="Coeff0" scale="1.0" value="+1.00000"   min="-10.0" max="10.0" free="0"/>
        <parameter name="Coeff1" scale="1.0" value="-0.1239176" min="-10.0" max="10.0" free="1"/>
        <parameter name="Coeff2" scale="1.0" value="+0.9751791" min="-10.0" max="10.0" free="1"/>
        <parameter name="Coeff3" scale="1.0" value="-3.0584577" min="-10.0" max="10.0" free="1"/>
        <parameter name="Coeff4" scale="1.0" value="+2.9089535" min="-10.0" max="10.0" free="1"/>
        <parameter name="Coeff5" scale="1.0" value="-1.3535372" min="-10.0" max="10.0" free="1"/>
        <parameter name="Coeff6" scale="1.0" value="+0.3413752" min="-10.0" max="10.0" free="1"/>
        <parameter name="Coeff7" scale="1.0" value="-0.0449642" min="-10.0" max="10.0" free="1"/>
        <parameter name="Coeff8" scale="1.0" value="+0.0024321" min="-10.0" max="10.0" free="1"/>
      </radialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>


CTA IRF background
^^^^^^^^^^^^^^^^^^

  .. code-block:: xml

    <source name="Background" type="CTAIrfBackground" instrument="CTA">
      <spectrum type="...">
        ...
      </spectrum>
    </source>


CTA cube background
^^^^^^^^^^^^^^^^^^^

  .. code-block:: xml

    <source name="Background" type="CTACubeBackground" instrument="CTA">
      <spatialModel type="GaussFunction">
        <parameter name="RA"    scale="1.0" value="83.6221" min="0.0"   max="360.0" free="0"/>
        <parameter name="DEC"   scale="1.0" value="22.01"   min="-90.0" max="90.0"  free="0"/>
        <parameter name="Sigma" scale="1.0" value="1.0"     min="0.0"   max="10.0"  free="0"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>

  .. code-block:: xml

    <source name="Background" type="CTACubeBackground" instrument="CTA">
      <spatialModel type="MapCubeFunction" file="map_cube.fits">
        <parameter name="Normalization" scale="1" value="1" min="0.001" max="1000.0" free="0"/>
      </spatialModel>
      <spectrum type="...">
        ...
      </spectrum>
    </source>



.. _sec_spectral_models:

Spectral model components
~~~~~~~~~~~~~~~~~~~~~~~~~

The following sections present the spectral model components that are available 
in ctools.

.. warning::

   The units of the spectral model depend on the spatial model component.
   Units quoted below apply to spatial models for a source component.
   Source intensities are generally given in units of
   :math:`{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}`.

   Exceptions to this rule exist for the following spatial models:

   - ``ConstantValue``: intensity units are given per steradian, e.g.
     :math:`{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}\,{\rm sr}^{-1}`

   - ``MapCubeFunction``: intensities are unitless, the spectral model 
     presents a relative scaling of the diffuse model

   If spectral models are used for a background model component, intensity 
   units are generally given in
   :math:`{\rm counts}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}\,{\rm sr}^{-1}`
   and correspond to the on-axis count rate.

   Exceptions to this rule exist for the following spatial models:

   - ``CTAIrfBackground``: intensities are unitless, the spectral model 
     presents a relative scaling of the background model

   - ``CTACubeBackground``: intensities are unitless, the spectral model 
     presents a relative scaling of the background model


Constant
^^^^^^^^

  .. code-block:: xml

   <spectrum type="ConstantValue">
     <parameter name="Normalization" scale="1e-16" value="5.7" min="1e-07" max="1000.0" free="1"/>
   </spectrum>

  The ``ConstantValue`` model component implements the constant function

  .. math::
    \frac{dN}{dE} = N_0

  where

  * :math:`N_0` = ``Normalization``
    :math:`[{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}]`



Power law
^^^^^^^^^

  .. code-block:: xml

   <spectrum type="PowerLaw">
     <parameter name="Prefactor" scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
     <parameter name="Index"     scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
     <parameter name="Scale"     scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
   </spectrum>

  The ``PowerLaw`` model component implements the power law function

  .. math::
    \frac{dN}{dE} = k_0 \left( \frac{E}{E_0} \right)^{\gamma}

  where

  * :math:`k_0` = ``Prefactor``
    :math:`[{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}]`
  * :math:`\gamma` = ``Index``
  * :math:`E_0` = ``Scale``
    :math:`[{\rm MeV}]`

  .. warning::

    The ``Scale`` parameter is not intended to be fitted.

  An alternative power law function that uses the integral flux as parameter
  rather than the Prefactor is specified by

  .. code-block:: xml

   <spectrum type="PowerLaw2">
     <parameter scale="1e-07" name="Integral"   min="1e-07" max="1000.0"    value="1.0" free="1"/>
     <parameter scale="1.0"   name="Index"      min="-5.0"  max="+5.0"      value="-2.0" free="1"/>
     <parameter scale="1.0"   name="LowerLimit" min="10.0"  max="1000000.0" value="100.0" free="0"/>
     <parameter scale="1.0"   name="UpperLimit" min="10.0"  max="1000000.0" value="500000.0" free="0"/>
   </spectrum>

  The ``PowerLaw2`` model component implements the power law function

  .. math::
    \frac{dN}{dE} = \frac{N(\gamma+1)E^{\gamma}}
                         {E_{\rm max}^{\gamma+1} - E_{\rm min}^{\gamma+1}}

  where

  * :math:`N` = ``Integral``
    :math:`[{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}]`
  * :math:`\gamma` = ``Index``
  * :math:`E_{\rm min}` = ``LowerLimit``
    :math:`[{\rm MeV}]`
  * :math:`E_{\rm max}` = ``UpperLimit``
    :math:`[{\rm MeV}]`

  .. warning::

    The ``LowerLimit`` and ``UpperLimit`` parameters are always treated as fixed and,
    the flux given by the ``Integral`` parameter is computed over the 
    range set by these two parameters.
    Use of this model allows the errors on the integral flux to be evaluated directly
    by :ref:`ctlike`.


Exponentially cut-off power law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  .. code-block:: xml

   <spectrum type="ExpCutoff">
     <parameter name="Prefactor" scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
     <parameter name="Index"     scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
     <parameter name="Cutoff"    scale="1e6"   value="1.0"  min="0.01"  max="1000.0" free="1"/>
     <parameter name="Scale"     scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
   </spectrum>

  The ``ExpCutoff`` model component implements the exponentially cut-off power law function

  .. math::
    \frac{dN}{dE} = k_0 \left( \frac{E}{E_0} \right)^{\gamma}
                    \exp \left( \frac{-E}{\tt E_{\rm cut}} \right)

  where

  * :math:`k_0` = ``Prefactor``
    :math:`[{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}]`
  * :math:`\gamma` = ``Index``
  * :math:`E_0` = ``Scale``
    :math:`[{\rm MeV}]`
  * :math:`E_{\rm cut}` = ``Cutoff``
    :math:`[{\rm MeV}]`

  .. warning::

    The ``Scale`` parameter is not intended to be fitted.



Broken power law
^^^^^^^^^^^^^^^^

  .. code-block:: xml

   <spectrum type="BrokenPowerLaw">
     <parameter name="Prefactor"  scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
     <parameter name="Index1"     scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
     <parameter name="BreakValue" scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="1"/>
     <parameter name="Index2"     scale="-1"    value="2.70" min="0.01"  max="1000.0" free="1"/>
   </spectrum>

  The ``BrokenPowerLaw`` model component implements the broken power law function

  .. math::

    \frac{dN}{dE} = k_0 \times \left \{
    \begin{eqnarray}
      \left( \frac{E}{E_b} \right)^{\gamma_1} & {\rm if\,\,} E < E_b \\
      \left( \frac{E}{E_b} \right)^{\gamma_2} & {\rm otherwise}
    \end{eqnarray}
    \right .

  where

  * :math:`k_0` = ``Prefactor``
    :math:`[{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}]`
  * :math:`\gamma_1` = ``Index1``
  * :math:`\gamma_2` = ``Index2``
  * :math:`E_b` = ``BreakValue``
    :math:`[{\rm MeV}]`

  .. warning::

    Note that the ``BreakValue`` parameter may be poorly constrained if 
    there is no clear spectral cut-off in the spectrum.
    This model may lead to complications in the maximum likelihood fitting.


Log parabola
^^^^^^^^^^^^

  .. code-block:: xml

   <spectrum type="LogParabola">
     <parameter name="Prefactor" scale="1e-17" value="5.878"   min="1e-07" max="1000.0" free="1"/>
     <parameter name="Index"     scale="-1"    value="2.32473" min="0.0"   max="+5.0"   free="1"/>
     <parameter name="Curvature" scale="-1"    value="0.074"   min="-5.0"  max="+5.0"   free="1"/>
     <parameter name="Scale"     scale="1e6"   value="1.0"     min="0.01"  max="1000.0" free="0"/>
   </spectrum>

  The ``LogParabola`` model component implements the log parabola function

  .. math::
    \frac{dN}{dE} = k_0 \left( \frac{E}{E_0} \right)^{\gamma+\eta \ln(E/E_0)}

  where

  * :math:`k_0` = ``Prefactor``
    :math:`[{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}]`
  * :math:`\gamma` = ``Index``
  * :math:`\eta` = ``Curvature``
  * :math:`E_0` = ``Scale``
    :math:`[{\rm MeV}]`

  .. warning::

    The ``Scale`` parameter is not intended to be fitted.

  An alternative XML format is supported for compatibility with the Fermi/LAT XML format:

  .. code-block:: xml

   <spectrum type="LogParabola">
     <parameter name="Prefactor" scale="1e-17" value="5.878"   min="1e-07" max="1000.0" free="1"/>
     <parameter name="alpha"     scale="1"     value="2.32473" min="0.0"   max="+5.0"   free="1"/>
     <parameter name="beta"      scale="1"     value="0.074"   min="-5.0"  max="+5.0"   free="1"/>
     <parameter name="Scale"     scale="1e6"   value="1.0"     min="0.01"  max="1000.0" free="0"/>
   </spectrum>

  where

  * ``alpha`` = -``Index``
  * ``beta`` = -``Curvature``


Gaussian
^^^^^^^^

  .. code-block:: xml

   <spectrum type="Gaussian">
     <parameter name="Normalization" scale="1e-10" value="1.0"  min="1e-07" max="1000.0" free="1"/>
     <parameter name="Mean"          scale="1e6"   value="5.0"  min="0.01"  max="100.0"  free="1"/>
     <parameter name="Sigma"         scale="1e6"   value="1.0"  min="0.01"  max="100.0"  free="1"/>
   </spectrum>

  The ``Gaussian`` model component implements the gaussian function

  .. math::
    \frac{dN}{dE} = \frac{N_0}{\sqrt{2\pi}\sigma}
                    \exp \left( \frac{-(E-\bar{E})^2}{2 \sigma^2} \right)

  where

  * :math:`N_0` = ``Normalization``
    :math:`[{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}]`
  * :math:`\bar{E}` = ``Mean``
    :math:`[{\rm MeV}]`
  * :math:`\sigma` = ``Sigma``
    :math:`[{\rm MeV}]`


File function
^^^^^^^^^^^^^

  .. code-block:: xml

   <spectrum type="FileFunction" file="data/filefunction.txt">
     <parameter scale="1.0" name="Normalization" min="0.0" max="1000.0" value="1.0" free="1"/>
   </spectrum>

  The ``FileFunction`` model component implements an arbitrary function 
  that is defined by intensity values at specific energies.
  The energy and intensity values are defined using an ASCII file with
  columns of energy and differential flux values.
  Energies are given in units of
  :math:`[{\rm MeV}]`,
  intensities are given in units of
  :math:`[{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}]`.
  The only parameter is a multiplicative normalization:

  .. math::
    \frac{dN}{dE} = N_0 \left. \frac{dN}{dE} \right\rvert_{\rm file}

  where

  * :math:`N_0` = ``Normalization``

  .. warning::

    If the file name is given as relative path, the path is relative to the
    working directory where the ctool has been launched.
    Alternatively, an absolute path may be specified.
    Any environment variable present in the path name will be expanded.


Node function
^^^^^^^^^^^^^

  .. code-block:: xml

   <spectrum type="NodeFunction">
     <node>
       <parameter name="Energy"    scale="1.0"   value="1.0" min="0.1"   max="1.0e20" free="0"/>
       <parameter name="Intensity" scale="1e-07" value="1.0" min="1e-07" max="1000.0" free="1"/>
     </node>
     <node>
       <parameter name="Energy"    scale="10.0"  value="1.0" min="0.1"   max="1.0e20" free="0"/>
       <parameter name="Intensity" scale="1e-08" value="1.0" min="1e-07" max="1000.0" free="1"/>
     </node>
   </spectrum>

  The ``NodeFunction`` model component implements a generalised broken 
  power law which is defined by a set of energy and intensity values
  (the so called nodes) that are piecewise connected by power laws.
  Energies are given in units of
  :math:`[{\rm MeV}]`,
  intensities are given in units of
  :math:`[{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}]`.

  .. warning::

    An arbitrary number of energy-intensity nodes can be defined in a node 
    function.
    The nodes need to be sorted by increasing energy.
    Although the fitting of the ``Energy`` parameters is formally possible 
    it may lead to numerical complications.
    If ``Energy`` parameters are to be fitted make sure that the ``min`` 
    and ``max`` attributes are set in a way that avoids inversion of the energy 
    ordering.


.. _sec_temporal_models:

Temporal model components
~~~~~~~~~~~~~~~~~~~~~~~~~

The following sections present the temporal model components that are available 
in ctools.

Constant
^^^^^^^^

  .. code-block:: xml

   <temporalModel type="Constant">
     <parameter name="Normalization" scale="1.0" value="1.0" min="0.1" max="10.0" free="0"/>
   </temporalModel>
