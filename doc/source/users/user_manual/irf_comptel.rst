.. _um_irf_comptel:

COMPTEL response functions
--------------------------

Formulation
~~~~~~~~~~~

A COMPTEL event is characterised by an instrument direction spanned by the
angles :math:`(\chi, \psi, \bar{\varphi})`. :math:`(\chi, \psi)` is the
direction of the photon after scattering in the upper detector plane, which is
determined from the photon interaction locations in both detector planes, and

.. math::
   \bar{\varphi} = \arccos \left( 1 - \frac{m_\mathrm{e}c^2}{E_2} + \frac{m_\mathrm{e}c^2}{E_1+E_2} \right)

is the Compton scattering angle as inferred from the energy deposits :math:`E_1`
and :math:`E_2` in the upper and lower detector planes, respectively.
The measured energy of the photon is estimated from the sum

.. math::
   E' = E_1 + E_2

of the energy deposits in both detector planes. The probability that a photon
which interacted in the upper detector plane will encounter a detector of the
lower plane is described by :math:`DRG(\chi, \psi, \bar{\varphi})`, which also
includes any cuts related to the removal of events coming from the Earth limb.

The COMPTEL response is factorised using

.. math::
   R(p',E',t'|p,E,t) = \frac{DRX(p)}{T} \times DRG(\chi, \psi, \bar{\varphi})
                       \times IAQ(\chi, \psi, \bar{\varphi} | p, E)

where
:math:`DRX(p)` is the exposure in units of :math:`cm^2 \, s`,
:math:`T` is the ontime in units of :math:`s`, and
:math:`IAQ(\chi, \psi, \bar{\varphi} | p, E)` quantifies the interaction
probability for a Compton scattering in the upper detector plane followed by
an interaction in the lower detector plane.


We note that :math:`IAQ(\chi, \psi, \bar{\varphi} | p, E)` is azimuthally
symmetric about the source direction, and the IAQ file is stored as a 2D FITS
image providing the interaction probabilities as function of
:math:`\bar{\varphi}` and :math:`\varphi_\mathrm{geo}` for a given energy
range, where :math:`\varphi_\mathrm{geo}` is the angular separation between
:math:`(\chi, \psi)` and :math:`p`.
