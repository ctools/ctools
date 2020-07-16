.. _um_models_bgd_iact:

IACT background models
----------------------

The following sections present the model components that are available in
ctools for the modelling of instrumental background in data from Imaging Air
Cherenkov Telescopes (IACTs) such as CTA, H.E.S.S., VERITAS, and MAGIC.

IACT background models are factorised into an optional spatial tag
(tags ``<spatialModel>`` or ``<radialModel>``) and a spectral tag
(tag ``<spectrum>``) using

.. math::
   M(p',E') = M_{\rm spatial}(p'|E') \times M_{\rm spectral}(E')

where :math:`M(p',E')` is given in units of
:math:`{\rm events} \,\, {\rm s}^{-1} {\rm MeV}^{-1} {\rm sr}^{-1}`.

For the spectral components, all spectra described in
:ref:`um_models_spectral`
may be used.

.. note::
   Spatial directions :math:`p'` and energies :math:`E'` are now the
   reconstructed quantities, hence no convolution of the model with the
   :ref:`Instrument Response Functions <um_irf>` is performed.


General IACT background
^^^^^^^^^^^^^^^^^^^^^^^

  The general IACT background model is factorised in a spatial and spectral
  component and has the type ``CTABackground``. It has the following XML
  structure:

  .. code-block:: xml

     <source name="Background" type="CTABackground" instrument="CTA">
       <spatialModel type="...">
       ...
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  .. warning::
     You need to specify the ``instrument`` label in the ``source`` tag
     that corresponds to the ``instrument`` label of the observation to which
     the model should apply. Supported ``instrument`` labels are ``CTA``,
     ``HESS``, ``VERITAS`` and ``MAGIC``.

     For example, if you analyse a H.E.S.S. observation you need to specify
     ``instrument="HESS"`` in the ``source`` tag, while the model type is
     still ``CTABackground``. So don't get confused!

     .. code-block:: xml

        <source name="Background" type="CTABackground" instrument="HESS">

     The same logic applies to the radial acceptance, IRF, effective area and
     cube background models.

  .. warning::
     In case that a background model should be used for the analysis of
     On/Off data, the corresponding ``OnOff`` ``instrument`` label needs to
     be selected. Supported ``OnOff`` ``instrument`` labels are ``CTAOnOff``,
     ``HESSOnOff``, ``VERITASOnOff`` and ``MAGICOnOff``.

     For example, if you analyse a H.E.S.S. On/Off observation you need to
     specify

     .. code-block:: xml

        <source name="Background" type="CTABackground" instrument="HESSOnOff">

     The same logic applies to the radial acceptance, IRF, effective area and
     cube background models.

  The following sections describe the spatial model components that are
  available.


Gaussian
~~~~~~~~

  The ``Gaussian`` model describes a 2D Gaussian shape in offset angle squared

  .. code-block:: xml

     <source name="Background" type="CTABackground" instrument="CTA">
       <spatialModel type="Gaussian">
         <parameter name="Sigma" scale="1.0" value="3.0" min="0.01" max="10.0" free="1"/>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  and implements

  .. math::
     M_{\rm spatial}(\theta) = \exp \left(-\frac{1}{2}
                               \left( \frac{\theta^2}{\sigma} \right)^2 \right)

  where

  * :math:`\sigma` = ``Sigma`` (degrees)

  and

  .. math::
     \theta = \sqrt{\mathrm{DETX} \times \mathrm{DETX} + \mathrm{DETY} \times\mathrm{DETY}}

  with :math:`\mathrm{DETX}` and :math:`\mathrm{DETY}` being the detector
  coordinates in the nominal system.


Profile
~~~~~~~

  The ``Profile`` model describes a radial profile

  .. code-block:: xml

     <source name="Background" type="CTABackground" instrument="CTA">
       <spatialModel type="Profile">
         <parameter name="Width" scale="1.0" value="1.5" min="0.1" max="1000.0" free="1"/>
         <parameter name="Core"  scale="1.0" value="3.0" min="0.1" max="1000.0" free="1"/>
         <parameter name="Tail"  scale="1.0" value="5.0" min="0.1" max="1000.0" free="1"/>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  and implements

  .. math::
     M_{\rm spatial}(\theta) = (1 + (\theta/c_0)^{c_1})^{-c_2/c_1}

  where

  * :math:`c_0` = ``Width`` (degrees)
  * :math:`c_1` = ``Core``
  * :math:`c_2` = ``Tail``


Polynom
~~~~~~~

  The ``Polynom`` model describes a polynomial with an arbitrary number of
  coefficients

  .. code-block:: xml

     <source name="Background" type="CTABackground" instrument="CTA">
       <spatialModel type="Polynom">
         <parameter name="Coeff0" scale="1.0" value="+1.00000"   min="-10.0" max="10.0" free="0"/>
         <parameter name="Coeff1" scale="1.0" value="-0.1239176" min="-10.0" max="10.0" free="1"/>
         <parameter name="Coeff2" scale="1.0" value="+0.9751791" min="-10.0" max="10.0" free="1"/>
         <parameter name="Coeff3" scale="1.0" value="-3.0584577" min="-10.0" max="10.0" free="1"/>
         ...
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  and implements

  .. math::
     M_{\rm spatial}(\theta) = \sum_{i=0}^m c_i \theta^i

  where

  * :math:`c_0` = ``Coeff0``
  * :math:`c_1` = ``Coeff1``
  * :math:`c_2` = ``Coeff2``
  * :math:`c_3` = ``Coeff3``
  * ...


Gradient
~~~~~~~~

  The ``Gradient`` model describes a bilinear gradient over the field of
  view

  .. code-block:: xml

     <source name="Background" type="CTABackground" instrument="CTA">
       <spatialModel type="Gradient">
         <parameter name="Grad_DETX" scale="1.0" value="0.0" min="-10.0" max="10.0" free="1"/>
         <parameter name="Grad_DETY" scale="1.0" value="0.0" min="-10.0" max="10.0" free="1"/>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  and implements

  .. math::
     M_{\rm spatial}(\mathrm{DETX},\mathrm{DETY}) =
     1 + \nabla_\mathrm{x} \mathrm{DETX} + \nabla_\mathrm{y} \mathrm{DETY}

  where

  * :math:`\nabla_\mathrm{x}` = ``Grad_DETX`` (per degree)
  * :math:`\nabla_\mathrm{y}` = ``Grad_DETY`` (per degree)



Multiplicative
~~~~~~~~~~~~~~

  The ``Multiplicative`` model describes a multiplication of spatial models

  .. code-block:: xml

     <source name="Background" type="CTABackground" instrument="CTA">
       <spatialModel type="Multiplicative">
         <spatialModel type="...">
           ...
         </spatialModel>
         <spatialModel type="...">
           ...
         </spatialModel>
         ...
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  and implements

  .. math::
     M_{\rm spatial}(\mathrm{DETX},\mathrm{DETY}) =
     \prod_{i=0}^{N-1} M^{(i)}_{\rm spatial}(\mathrm{DETX},\mathrm{DETY})

  where :math:`M^{(i)}_{\rm spatial}(\mathrm{DETX},\mathrm{DETY})` is any
  spatial model component, including another multiplicative model, and
  :math:`N` is the number of model components that are multiplied.
  For example, the default model for a H.E.S.S. data analysis is specified
  by

  .. code-block:: xml

     <source name="Background" type="CTABackground" instrument="CTA">
       <spatialModel type="Multiplicative">
         <spatialModel type="Gaussian">
           <parameter name="Sigma" scale="1.0" value="3.0" min="0.01" max="10.0" free="1"/>
         </spatialModel>
         <spatialModel type="Gradient">
           <parameter name="Grad_DETX" scale="1.0" value="0.0" min="-10.0" max="10.0" free="1"/>
           <parameter name="Grad_DETY" scale="1.0" value="0.0" min="-10.0" max="10.0" free="1"/>
         </spatialModel>
       </spatialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>


Radial acceptance background
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  For legacy reasons, there exists a class of radially symmetric background
  models of the type ``RadialAcceptance`` with the following XML structure:

  .. code-block:: xml

     <source name="Background" type="RadialAcceptance" instrument="CTA">
       <radialModel type="Gaussian">
         ...
       </radialModel>
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  These models require a ``<radialModel>`` tag as the spatial component and
  accept all spatial model types that take the offset angle :math:`\theta`
  as variable, such as ``Gaussian``, ``Profile`` and ``Polynom``.

  .. warning::
     The use of the radial acceptance model is deprecated, and the
     ``CTABackground`` model should be used instead.


IRF background
^^^^^^^^^^^^^^

  The :ref:`Instrument Response Functions (IRFs) <um_irf>` contain a template
  that predicts the background rate as function of position in the field of view
  and measured energy in units of
  :math:`{\rm events} \, {\rm s}^{-1} {\rm MeV}^{-1} {\rm sr}^{-1}`. This template
  can be used by specifying a model of type ``CTAIrfBackground``. No spatial component
  will be specified explicitly since the spatial (and spectral) information is
  already contained in the template.

  The model will be multiplied by a spectral component to allow for the adjustment
  of the energy distribution of the background rate.

  .. code-block:: xml

     <source name="Background" type="CTAIrfBackground" instrument="CTA">
       <spectrum type="...">
         ...
       </spectrum>
     </source>

  If the observation is an On/Off observation, do not forget to switch the instrument
  to ``CTAOnOff``:

  .. code-block:: xml

     <source name="Background" type="CTAIrfBackground" instrument="CTAOnOff">
       <spectrum type="...">
         ...
       </spectrum>
     </source>


Effective area background
^^^^^^^^^^^^^^^^^^^^^^^^^

  Instead of using the background template the effective area for gamma rays can
  also be used to model the instrumental background. Note that in this case the
  effective area has to be scaled to a reasonable background rate by selecting
  appropriate values for the spectral model component.

  .. code-block:: xml

     <source name="Background" type="CTAAeffBackground" instrument="CTA">
       <spectrum type="...">
         ...
       </spectrum>
     </source>


Cube background
^^^^^^^^^^^^^^^

  For a stacked analysis, the background rates are predicted by a so called
  background cube. The FITS file name of the background cube is specified
  either as ctools task parameter, or using the ``BkgCube`` parameter in the
  the observation definition XML file.

  The background cube model is used by specifying a model of type
  ``CTACubeBackground``. Similar to the ``CTAIrfBackground`` model, the
  background cube is multplied with a spectral model to allow for the adjustment
  of the energy distribution of the background rate.

  .. code-block:: xml

     <source name="Background" type="CTACubeBackground" instrument="CTA">
       <spectrum type="...">
         ...
       </spectrum>
     </source>
