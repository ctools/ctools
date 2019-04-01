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


Event types
~~~~~~~~~~~

The LAT events are partitioned into exclusive event types that for Pass 6 and
Pass 7 data correspond to pair conversions located in either the front or the
back section of the tracker. For Pass 8 the event partitioning has been
generalised to other event types. For each event type, a specific response
function exists that will be designated in the following with the superscript
:math:`\alpha`.


Livetime cube
~~~~~~~~~~~~~

The livetime cube is a means to speed up the exposure calculations in a
Fermi/LAT analysis and contains the integrated livetime as a function of sky
position and inclination angle with respect to the LAT z-axis.
This livetime, denoted by :math:`\tau(p,\theta)`, is the time that the LAT
observed a given position on the sky at a given inclination angle, and includes
the history of the LAT's orientation during the entire observation period.
A Fermi/LAT livetime cube includes also a version of the livetime information
that is weighted by the livetime fraction (i.e. the ratio between livetime and
ontime) and that allows correction of inefficiencies introduced by so-called
ghost events, and that we denote here by :math:`\tau_\mathrm{wgt}(p,\theta)`.


Mean point-source PSF
~~~~~~~~~~~~~~~~~~~~~

GammaLib, the library that is underlying ctools, natively implements the
computation of the mean PSF for point sources.
The exposure for a given sky direction :math:`p`, photon energy :math:`E` and
event type :math:`\alpha` is computed using

.. math::
   X^\alpha(p, E) = f_1^\alpha(E) \int_{\theta} \tau(p,\theta) \,
                    A_\mathrm{eff}^\alpha(E,\theta) \, d\theta
                  + f_2^\alpha(E) \int_{\theta} \tau_\mathrm{wgt}(p,\theta)
                    \, A_\mathrm{eff}^\alpha(E,\theta) \, d\theta

The exposure weighted point spread function is computed using

.. math::
   \mathrm{\it PSF}^\alpha(\delta|p,E) =
          f_1^\alpha(E) \int_{\theta} \tau(p,\theta) \,
          A_\mathrm{eff}^\alpha(E, \theta) \,
          \mathrm{\it PSF}^\alpha(\delta|E,\theta) \, d\theta
        + f_2^\alpha(E) \int_{\theta} \tau_\mathrm{wgt}(p,\theta) \,
          A_\mathrm{eff}^\alpha(E,\theta) \, \mathrm{\it PSF}^\alpha(\delta|E,\theta) \,
          d\theta,

where :math:`f_1^\alpha(E)` and :math:`f_2^\alpha(E)` are energy and event
type dependent efficiency factors.

The mean point spread function for a point source is computed using

.. math::
   \overline{\mathrm{\it PSF}}(\delta|p,E) =
             \frac{\sum_\alpha \mathrm{\it PSF}^\alpha(\delta|p,E)}
                  {\sum_\alpha X^\alpha(p,E)}

where the sum is taken over all event types :math:`\alpha`.

