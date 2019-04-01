.. _um_classical:

Classical On/Off analysis
-------------------------

The traditional technique that is widely used by the IACT community for spectral
analysis consists in selecting the source events from an On region and
estimating the background events from one or several signal-free Off regions.
The events are put into On and Off spectra, and the effective response is
computed for the On region that allows to turn the event spectrum into flux
points.

The cscript :ref:`csphagen` generates all files that are necessary for an
On/Off spectral analysis and saves them in the
`OGIP <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node5.html>`_
format that is normally used in X-ray astronomy and that is compliant with the
`XSPEC <https://heasarc.gsfc.nasa.gov/xanadu/xspec/>`_
spectral fitting package. This format is composed of Pulse Height Analyzer
spectral files (``PHA``), an Auxiliary Response File (``ARF``) and a Redistribution
Matrix File (``RMF``).

``PHA`` files are generated for the On and the Off region by binning the events
in both regions as function of reconstructed energy and storing the corresponding
vectors :math:`n^\mathrm{on}_i` and :math:`n^\mathrm{off}_i`, where :math:`i`
is the index of the energy bin.

The ``ARF`` is computed using

.. math::
   ARF(E) = \int_\mathrm{on} \int_{p} A_{\rm eff}(p,E) \times PSF(p'|p,E)
            \times M_s(p,E) \, dp \, dp'

where
:math:`A_{\rm eff}(p,E)` is the effective area,
:math:`PSF(p'|p,E)` is the point spread function,
:math:`M_s(p,E)` is the source model,
:math:`p` and :math:`p'` are the true and reconstructed photon arrival
directions, respectively, and :math:`E` is the true energy. The integration
in :math:`p'` is done over the On region, the integration in :math:`p` is done
over all :math:`p` that contribute events within the On region. These
integrations assure the correct computation of the number of source events
within the On region, even if the emission is not fully contained within that
region.

.. note::
   The source model is normalised to unity,

   .. math::
      \int_{p} M_s(p,E) \, dp = 1

   where the integral over :math:`p` is taken over the full sky.

.. note::
   The convolution with the PSF is skipped and the model :math:`M_s(p,E)`
   ignored if the IRFs are calculated for the events surviving a directional
   cut around the assumed source position, which is specified by the ``RAD_MAX``
   keyword in the IRF files. The user must take care of specifying an On region
   compatible with this directional cut. In that case

   .. math::
      ARF(E) = \frac{\int_\mathrm{on} A_{\rm eff}(p',E) \, dp'}
                    {\int_\mathrm{on} dp'}

   which is the mean effective area over the On region.

The ``RMF`` is computed using

.. math::
   RMF_i(E) = \frac{\int_\mathrm{on} \int_{E'} A_{\rm eff}(p,E) \times
                    E_{\rm disp}(E'|p,E) \, dE' \, dp}
                    {\int_\mathrm{on} A_{\rm eff}(p,E) \, dp}

where :math:`E_{\rm disp}(E'|p,E)` is the energy dispersion. The integration in
:math:`p` is done over the On region to accommodate for possible variations of
the energy dispersion over that region. The integration over reconstructed
energy :math:`E'` is done over the width of the energy bin :math:`i`.

:ref:`csphagen` also computes the background scaling factors :math:`\alpha_i`
and background response vectors :math:`b_i`.

The background scaling factors :math:`\alpha_i` are stored in the
``BACKSCAL`` column of the On PHA file and are computed using

.. math::
   \alpha_i = \frac{\int_\mathrm{on} M_b(p',E') \, dp'}
                   {\int_\mathrm{off} M_b(p',E') \, dp'}

where :math:`M_b(p',E')` is a background acceptance model, specified either
using a model definition XML file, or the template background found in the IRF.
If no background acceptance model is provided, :math:`M_b(p',E')=1`, and
:math:`\alpha_i` gives the solid angle ratio between On and Off regions.

The background response vectors :math:`b_i` are stored in the ``BACKRESP``
column of each Off PHA file and are computed using

.. math::
   b_i = \int_\mathrm{off} M_b(p',E') \, dp'

where :math:`M_b(p',E')` are evaluated at the reconstructed energy bins
:math:`i`. If no background acceptance model is provided, :math:`M_b(p',E')=1`,
and :math:`b_i` gives the solid angle of the Off region.



