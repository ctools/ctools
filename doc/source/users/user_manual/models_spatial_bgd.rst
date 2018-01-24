.. _sec_models_spatial_bgd:

Spatial background model components
-----------------------------------

The following sections present the spatial model components that are available 
in ctools for instrumental background modelling for CTA.

CTA radial background
^^^^^^^^^^^^^^^^^^^^^

  There exist a number of radial CTA background models that factorise into a
  spatial and spectral component, and that model the spatial component as a
  function that only dependends on the offset angle :math:`\theta` from the
  pointing direction. These models have the type ``RadialAcceptance`` and
  require a ``<radialModel>`` tag as the spatial component.

  The first radial model is of type ``Gaussian``

  .. code-block:: xml

     <source name="Background" type="RadialAcceptance" instrument="CTA">
       <radialModel type="Gaussian">
         <parameter name="Sigma" scale="1.0" value="3.0" min="0.01" max="10.0" free="1"/>
       </radialModel>
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1e-6" value="61.8" min="0.0"  max="1000.0" free="1"/>
         <parameter name="Index"       scale="-1"   value="1.85" min="0.0"  max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
       </spectrum>
     </source>

  and implements

  .. math::
     M_{\rm spatial}(\theta) = \exp \left(-\frac{1}{2}
                               \left( \frac{\theta^2}{\sigma} \right)^2 \right)

  where

  * :math:`\sigma` = ``Sigma`` (degrees)

  The second radial model is of type ``Profile``

  .. code-block:: xml

     <source name="Background" type="RadialAcceptance" instrument="CTA">
       <radialModel type="Profile">
         <parameter name="Width" scale="1.0" value="1.5" min="0.1" max="1000.0" free="1"/>
         <parameter name="Core"  scale="1.0" value="3.0" min="0.1" max="1000.0" free="1"/>
         <parameter name="Tail"  scale="1.0" value="5.0" min="0.1" max="1000.0" free="1"/>
       </radialModel>
        <spectrum type="PowerLaw">
          <parameter name="Prefactor"   scale="1e-6" value="61.8" min="0.0"  max="1000.0" free="1"/>
          <parameter name="Index"       scale="-1"   value="1.85" min="0.0"  max="+5.0"   free="1"/>
          <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
        </spectrum>
     </source>

  and implements

  .. math::
     M_{\rm spatial}(\theta) = (1 + (\theta/c_0)^{c_1})^{-c_2/c_1}

  where

  * :math:`c_0` = ``Width`` (degrees)
  * :math:`c_1` = ``Core``
  * :math:`c_2` = ``Tail``

  The third radial model is of type ``Polynom``

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
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1e-6" value="61.8" min="0.0"  max="1000.0" free="1"/>
         <parameter name="Index"       scale="-1"   value="1.85" min="0.0"  max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
       </spectrum>
     </source>

  and implements

  .. math::
     M_{\rm spatial}(\theta) = \sum_{i=0}^m c_i \theta^i

  where

  * :math:`c_0` = ``Coeff0``
  * :math:`c_1` = ``Coeff1``
  * ...

  (the number of polynomial coefficients is arbitrary).


CTA IRF background
^^^^^^^^^^^^^^^^^^

  The :ref:`Instrument Response Functions (IRFs) <sec_response>` contain a template
  that predicts the background rate as function of position in the field of view
  and measured energy in units of
  :math:`{\rm events} \, {\rm s}^{-1} {\rm MeV}^{-1} {\rm sr}^{-1}`. This template
  can be used by specifying a model of type ``CTAIrfBackground``. No spatial component
  will be specified explicitly since the spatial (and spectral) information is
  already contained in the template. The model will be multiplied by a spectral law.
  In the example below, the template is multiplied by a power law with normalization
  of 1 and slope 0 (i.e. the template is taken as is, but a fit can adjust the model
  to compensate for inaccuracies):

  .. code-block:: xml

     <source name="Background" type="CTAIrfBackground" instrument="CTA">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1.0" value="1.0" min="1e-3" max="1e3"    free="1"/>
         <parameter name="Index"       scale="1.0" value="0.0" min="-5.0" max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6" value="1.0" min="0.01" max="1000.0" free="0"/>
       </spectrum>
     </source>

  If the observation is an On/Off observation, do not forget to switch the instrument
  to ``CTAOnOff``:

  .. code-block:: xml

     <source name="Background" type="CTAIrfBackground" instrument="CTAOnOff">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1.0" value="1.0" min="1e-3" max="1e3"    free="1"/>
         <parameter name="Index"       scale="1.0" value="0.0" min="-5.0" max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6" value="1.0" min="0.01" max="1000.0" free="0"/>
       </spectrum>
     </source>


CTA effective area background
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  Instead of using the background template the effective area for gamma rays can
  also be used to model the instrumental background. Note that in this case the
  effective area has to be scaled to a reasonable background rate. An example is
  given below:

  .. code-block:: xml

     <source name="Background" type="CTAAeffBackground" instrument="CTA">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1e-14" value="1.0"  min="1e-3" max="1e3"    free="1"/>
         <parameter name="Index"       scale="1.0"   value="-2.4" min="-5.0" max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6"   value="1.0"  min="0.01" max="1000.0" free="0"/>
       </spectrum>
     </source>


CTA cube background
^^^^^^^^^^^^^^^^^^^

  For a stacked analysis, the background rate is predicted by a so called background
  cube. The background cube is used by specifying a model of type ``CTACubeBackground``.
  Similar to the ``CTAIrfBackground`` model, the background cube is multplied by
  a spectral model.

  .. code-block:: xml

     <source name="Background" type="CTACubeBackground" instrument="CTA">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1.0" value="1.0" min="1e-3" max="1e3"    free="1"/>
         <parameter name="Index"       scale="1.0" value="0.0" min="-5.0" max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6" value="1.0" min="0.01" max="1000.0" free="0"/>
       </spectrum>
     </source>
