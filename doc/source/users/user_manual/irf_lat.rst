.. _um_irf_lat:

Fermi/LAT response functions
----------------------------

Formulation
~~~~~~~~~~~

The instrument response functions for Fermi/LAT are factorised into
the effective area :math:`A_{\rm eff}(p,\theta)` (units :math:`cm^2`),
the point spread function :math:`PSF(\delta|E,\theta)`,
and the energy dispersion :math:`E_{\rm disp}(E'|E,\theta)`
following

.. math::
    R(p',E',t'|p,E,t) = A_{\rm eff}(p,\theta) \times PSF(\delta|E,\theta)
                        \times E_{\rm disp}(E'|E,\theta)

where :math:`\theta` is the inclination angle with respect to the LAT z-axis,
and :math:`\delta` is the angular separation between the true and measured
photon directions :math:`p` and :math:`p'`, respectively,

.. math::
   \int PSF(p'|p,E,t) \, dp' = 1

and

.. math::
   \int E_{\rm disp}(E'|p,E,t) \, dE' = 1

The instrument response function is independent of time.

Two functional forms are available for the point spread function which are both
composed of a superposition of two King functions:

.. math::
   \mathrm{\it PSF}_1(\delta | E, \theta) =
   n_\mathrm{c}
   \left( 1-\frac{1}{\gamma_\mathrm{c}} \right)
   \left( 1 + \frac{1}{2\gamma_\mathrm{c}} \frac{\delta^2}{\sigma^2} \right)^{-\gamma_\mathrm{c}}
   + n_\mathrm{t}
   \left( 1-\frac{1}{\gamma_\mathrm{t}} \right)
   \left( 1 + \frac{1}{2\gamma_\mathrm{t}} \frac{\delta^2}{\sigma^2} \right)^{-\gamma_\mathrm{t}}

and

.. math::
   \mathrm{\it PSF}_3(\delta | E, \theta) =
   n_\mathrm{c}
   \left(
   \left( 1-\frac{1}{\gamma_\mathrm{c}} \right)
   \left( 1 + \frac{1}{2\gamma_\mathrm{c}} \frac{\delta^2}{s_\mathrm{c}^2} \right)^{-\gamma_\mathrm{c}} \right.
   \left.
   + n_\mathrm{t}
   \left( 1-\frac{1}{\gamma_\mathrm{t}} \right)
   \left( 1 + \frac{1}{2\gamma_\mathrm{t}} \frac{\delta^2}{s_\mathrm{t}^2} \right)^{-\gamma_\mathrm{t}}
   \right).

The parameters :math:`n_\mathrm{c}`, :math:`n_\mathrm{t}`, :math:`s_\mathrm{c}`,
:math:`s_\mathrm{t}`, :math:`\sigma`, :math:`\gamma_\mathrm{c}`,
:math:`\gamma_\mathrm{t}` depend on energy :math:`E` and off-axis angle
:math:`\theta`.
Energy dispersion is so far not implemented.
