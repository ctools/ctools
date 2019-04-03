.. _um_models_temporal:

Temporal model components
-------------------------

The following sections present the temporal model components that are available 
in ctools.

.. warning::
   Recall that the model is factorised according to

   .. math::
      M(p,E,t) = M_{\rm spatial}(p|E) \times M_{\rm spectral}(E) \times M_{\rm temporal}(t)

   which implies that the temporal model component is multiplied with the
   spectral model component.


Constant
^^^^^^^^

  .. code-block:: xml

     <temporal type="Constant">
       <parameter name="Normalization" scale="1.0" value="1.0" min="0.1" max="10.0" free="0"/>
     </temporal>

  This temporal model component implements a constant source

  .. math::
     M_{\rm temporal}(t) = N_0

  where

  * :math:`N_0` = ``Normalization``


Light Curve
^^^^^^^^^^^

  .. code-block:: xml

     <temporal type="LightCurve" file="model_temporal_lightcurve.fits">
       <parameter name="Normalization" scale="1" value="1.0" min="0.0" max="1000.0" free="0"/>
     </temporal>

  This temporal model component implements a light curve :math:`r(t)`

  .. math::
     M_{\rm temporal}(t) = N_0 \times r(t)

  where

  * :math:`N_0` = ``Normalization``

  The light curve is defined by nodes in a FITS file that specify the relative
  flux normalization as function of time (file ``model_temporal_lightcurve.fits``
  in the example above). The structure of the light curve FITS
  file is shown in the figure below. The light curve is defined in the first
  extension of the FITS file and consists of a binary table with the columns
  ``TIME`` and ``NORM``. Times in the ``TIME`` columns are given in seconds
  and are counted with respect to a time reference that is defined in the
  header of the binary table. Times need to be specified in ascending order.
  The values in the ``NORM`` column specify :math:`r(t)` at times :math:`t`,
  and should be comprised between 0 and 1.

  .. _fig_model_lightcurve:

  .. figure:: models_lightcurve.png
     :align: center
     :width: 100%

     *Structure of light curve FITS file*

  .. warning::
     Fitting of light curves only makes sense for an unbinned maximum likelihood
     analysis, since in a binned or stacked analysis the times of individual
     events are dropped.


Phase Curve
^^^^^^^^^^^

  .. code-block:: xml

     <temporal type="PhaseCurve" file="model_temporal_phasecurve.fits">
       <parameter name="Normalization" scale="1" value="1.0"     min="0.0" max="1000.0"   free="0"/>
       <parameter name="MJD"           scale="1" value="51544.5" min="0.0" max="100000.0" free="0"/>
       <parameter name="Phase"         scale="1" value="0.0"     min="0.0" max="1.0"      free="0"/>
       <parameter name="F0"            scale="1" value="1.0"     min="0.0" max="1000.0"   free="0"/>
       <parameter name="F1"            scale="1" value="0.1"     min="0.0" max="1000.0"   free="0"/>
       <parameter name="F2"            scale="1" value="0.01"    min="0.0" max="1000.0"   free="0"/>
     </temporal>

  This temporal model component implements a phase curve :math:`r(\Phi(t))`

  .. math::
     M_{\rm temporal}(t) = N_0 \times r(\Phi(t))

  where the phase as function of time is computed using

  .. math::
     \Phi(t) = \Phi_0 + f(t-t_0) + \frac{1}{2}\dot{f} (t-t_0)^2 +
                                   \frac{1}{6}\ddot{f} (t-t_0)^3

  and

  * :math:`N_0` = ``Normalization``
  * :math:`t_0` = ``MJD``
  * :math:`\Phi_0` = ``Phase``
  * :math:`f` = ``F0``
  * :math:`\dot{f}` = ``F1``
  * :math:`\ddot{f}` = ``F2``

  The phase curve is defined by nodes in a FITS file that specify the relative
  flux normalization as function of phase (file ``model_temporal_phasecurve.fits``
  in the example above). The structure of the phase curve
  FITS file is shown in the figure below. The phase curve is defined in the
  first extension of the FITS file and consists of a binary table with the
  columns ``PHASE`` and ``NORM``. Phase values in the ``PHASE`` column need to
  be comprised between 0 and 1 and need to be given in ascending order. The
  values in the ``NORM`` column specify :math:`r(\Phi(t))` at phases
  :math:`\Phi(t)`, and should be comprised between 0 and 1.

  .. _fig_models_phasecurve:

  .. figure:: models_phasecurve.png
     :align: center
     :width: 40%

     *Structure of phase curve FITS file*

  By default, the ``NORM`` values are recomputed internally so that the
  phase-averaged normalisation is one, i.e.

  .. math::
     \int_0^1 r(\Phi) d\Phi = 1

  In that case, the spectral component corresponds to the phase-averaged
  spectrum. If the internal normalisation should be disabled the
  ``normalize="0"`` attribute needs to be added to the temporal tag, i.e.

  .. code-block:: xml

     <temporal type="PhaseCurve" file="model_temporal_phasecurve.fits" normalize="0">

  In that case the ``NORM`` values are directly multiplied with the spectral
  component.

  .. warning::
     Fitting of phase curves only makes sense for an unbinned maximum likelihood
     analysis, since in a binned or stacked analysis the times of individual
     events are dropped.

  .. warning::
     Fitting of phase curve parameters may not properly work for pulsar
     frequencies.
