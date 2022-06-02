.. _um_models_source:

Source models
-------------

The general model is describe in ctools using a
:ref:`model definition XML file <glossary_moddef>`.
Below is a simple example of such a file comprising one source and one 
background model (see :ref:`um_models_bgd_iact`, :ref:`um_models_bgd_comptel`
and :ref:`um_models_bgd_spi` for the implementation of the
background models).

Each source model is factorised into a spatial (tags ``<spatialModel>``, a
spectral (tag ``<spectrum>``), and a temporal component (tag ``<temporal>``)
using

.. math::
   M(p,E,t) = M_{\rm spatial}(p|E) \times M_{\rm spectral}(E) \times M_{\rm temporal}(t)

where :math:`M(p,E,t)` is given in units of
:math:`photons \,\, {\rm cm}^{-1} {\rm s}^{-1} {\rm MeV}^{-1} {\rm sr}^{-1}`.

The spatial model component describes the energy dependent morphology of the
source.
It satisfies

.. math::
   \int_{\Omega} M_{\rm spatial}(p|E) \, d\Omega = 1

for all :math:`E`, hence the spatial component does not
impact the spatially integrated spectral and temporal properties of the
source (the integration is done here over the spatial parameters
:math:`p` in a spherical coordinate system).
The units of the spatial model component are
:math:`[M_{\rm spatial}] = {\rm sr}^{-1}`
(but see the exception of the
:ref:`DiffuseMapCube spatial model component <um_models_spatial_src_diffuse_cube>`).

The spectral model component describes the spatially integrated spectral
distribution of the source.
It satisfies

.. math::
   \int_{E} M_{\rm spectral}(E) \, dE = \Phi

where :math:`\Phi` is the spatially and spectrally integrated total source
flux.
The units of the spectral model component
are :math:`[M_{\rm spectral}] = {\rm cm}^{-2} {\rm s}^{-1} {\rm MeV}^{-1}`.

The temporal model component describes the relative variation of the
source flux with respect to the mean value given by the spectral model
component.
The temporal model component is unit less.

In this specific example, the source component ``Crab`` describes 
a point source at the location of the Crab nebula with a power law spectral
shape that is constant in time.
The background component ``Background`` is modelled using the template that is
provided with the
:ref:`Instrument Response Functions (IRFs) <um_irf>`
and that is multiplied by a spectral power law function.

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <source_library title="source library">
    <source name="Crab" type="PointSource">
      <spatialModel type="PointSource">
        <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
        <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
      </spatialModel>
      <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"       scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
      </spectrum>
      <temporal type="Constant">
        <parameter name="Normalization" scale="1.0" value="1.0" min="0.0" max="1000.0" free="0"/>
      </temporal>
    </source>
    <source name="CTABackgroundModel" type="CTAIrfBackground" instrument="CTA">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor"   scale="1.0"  value="1.0"  min="1e-3" max="1e+3"   free="1"/>
        <parameter name="Index"       scale="1.0"  value="0.0"  min="-5.0" max="+5.0"   free="1"/>
        <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
      </spectrum>
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
   For compatibility reasons with the Fermi/LAT XML format the ``<temporal>``
   tag can be omitted for models that are constant in time.

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

.. note::

   The syntax of the :ref:`model definition XML file <glossary_moddef>` has been
   inspired from the syntax used by the Fermi/LAT ScienceTools, but for reasons
   of clarity and homogenity of the various model and parameter names we have
   made some modifications.
   Nevertheless, the format used by the Fermi/LAT ScienceTools is also
   supported.

.. note::

   To cope with cross-calibration uncertainties between different instruments,
   instrument specific scaling factors :math:`S_{\rm inst}` can be multiplied
   to the source models according to

   .. math::
      M(p,E,t) = S_{\rm inst} \times M_{\rm spatial}(p|E) \times M_{\rm spectral}(E) \times M_{\rm temporal}(t)

   and an example for the specification of scale factors :math:`S_{\rm CTA}` for
   ``CTA`` and :math:`S_{\rm HESS}` for H.E.S.S. in the
   :ref:`model definition XML file <glossary_moddef>`
   is shown below.
   In the example, :math:`S_{\rm CTA}` is fixed to :math:`1`, hence the
   scale factor could also be omitted from the XML file.
   The scale factor :math:`S_{\rm HESS}` is set to :math:`1.1` and will be
   fitted during the model fit.

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
          <scaling>
            <instrument name="CTA"  scale="1.0" min="0.1" max="10.0" value="1.0" free="0"/>
            <instrument name="HESS" scale="1.0" min="0.1" max="10.0" value="1.1" free="1"/>
          </scaling>
        </source>
      </source_library>

