.. _um_likelihood:

Maximum likelihood estimation
-----------------------------

The central method behind the ctools data analysis is the maximum likelihood
estimation of the parameters of a given :ref:`model <um_models>`. The method
obtains the parameter estimates by finding the parameter values that maximize
the likelihood function :math:`L(M)`. The likelihood function quantifies the
probability that the data are drawn from a particular :ref:`model <um_models>`
:math:`M`. The formulae used for the likelihood computation depends on whether
the data are unbinned or binned and on the assumed underlying statistical law.
An iterative Levenberg-Marquardt algorithm is employed for maximum likelihood
estimation.


Unbinned data
~~~~~~~~~~~~~

For unbinned data the Poisson formula

.. math::
   -\ln L(M) = E(M) - \sum_i \ln P(p'_i, E'_i ,t'_i | M)

is used, where the sum is taken over all events :math:`i`, characterised by the
instrument direction :math:`p'_i`, the measured energy :math:`E'_i`
and the trigger time :math:`t'_i`. :math:`P(p', E' ,t' | M)` is
the probability density that given the :ref:`model <um_models>` :math:`M`, an
event with instrument direction :math:`p'`, measured energy :math:`E'`
and trigger time :math:`t'` occurs. :math:`E(M)` is the total number of events
that are predicted to occur during an observation given the
:ref:`model <um_models>` :math:`M`, computed by integrating the probability
density over the trigger time, measured energy and instrument direction:

.. math::
   E(M) = \int_{GTI} \int_{Ebounds} \int_{ROI} P(p',E',t' | M) \,
   d\vec{p}' \, dE' \, dt'

The temporal integration boundaries are defined by the Good Time Intervals
(GTIs) that define contiguous periods in time during which data were taken.
The spatial integration boundaries are defined by a so-called Region of
Interest (ROI).


Binned data
~~~~~~~~~~~

For binned data following a Poisson distribution the formula

.. math::
   -\ln L(M) = \sum_i e_i(M) - n_i \ln e_i(M)

is used, where the sum over :math:`i` is now taken over all data cube bins.
:math:`n_i` is the number of events observed in bin :math:`i`, and

.. math::
   e_i(M) = P(p'_i, E'_i ,t'_i | M) \times \Omega_i \times \Delta E_i \times
   \Delta T_i

is the predicted number of events from :ref:`model <um_models>` :math:`M` in
bin :math:`i`. The probability density is evaluated for the reference
instrument direction :math:`p'_i`, measured energy :math:`E'_i` and trigger
time :math:`t'_i` of bin :math:`i`, taken to be the values at the bin centre,
and multiplied by the solid angle :math:`\Omega_i`, the energy width
:math:`\Delta E_i` and the exposure time (or ontime) :math:`\Delta T_i` of
bin :math:`i` (note that a live time correction is included in the probability
density :math:`P(p'_i, E'_i ,t'_i | M)` so that the multiplication is done by
the true exposure time).


Binned data following Gaussian statistic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, if the data follow a Gaussian distribution the formula

.. math::
   -\ln L(M) = \frac{1}{2} \sum_i \left( \frac{n_i - e_i(M)}{\sigma_i} \right)^2

is used, where :math:`\sigma_i` is the statistical uncertainty in the measured
number of events for bin :math:`i`.


On/Off analysis
~~~~~~~~~~~~~~~

For a so-called On/Off analysis, where spectra from a source (or On) region
and source-free (or Off) region are analysed jointly, there exist two options:
``cstat`` and ``wstat``.

cstat
^^^^^

``cstat`` makes use of a background model :math:`M_b` and employs the Poisson
formula

.. math::
   -\ln L(M) = \sum_i s_i(M_s) + \alpha_i(M_b) \, b_i(M_b) -
               n^\mathrm{on}_i \ln [s_i(M_s)+ \alpha_i(M_b) \, b_i(M_b)] + b_i(M_b)  -
               n^\mathrm{off}_i \ln b_i(M_b)

where the :ref:`model <um_models>` :math:`M` is split into signal (i.e. gamma
rays) and background, i.e. :math:`M = M_s + M_b`,
:math:`n^\mathrm{on}_i` is the number of events in bin :math:`i` of the On
region,
:math:`n^\mathrm{off}_i` is the number of events in bin :math:`i` of the Off
region,
:math:`s_i(M_s)` is the number of expected signal counts in bin :math:`i` of
the On region,
:math:`b_i(M_b)` is the number of expected background counts in bin :math:`i`
of the Off region,
and

.. math::
   \alpha_i(M_b) = \frac{\int_\mathrm{on} M_b d\Omega}{\int_\mathrm{off} M_b d\Omega}

is the ratio between the spatial integral over the background model in the On
region and the Off region for bin :math:`i`.

wstat
^^^^^

``wstat`` does not make use of an explicit background model but assumes that
the background rate per solid angle is the same in the On and the Off region.
The Poisson formula is then

.. math::
   -\ln L(M_s) = \sum_i s_i(M_s) + \alpha_i b_i(M_s) -
                 n^\mathrm{on}_i \ln  [s_i(M_s)+ \alpha_i \, b_i(M_s)] + b_i(M_s) -
                 n^\mathrm{off}_i \ln b_i(M_s) -\\
                 n^\mathrm{on}_i (1-\ln n^\mathrm{on}_i) -
                 n^\mathrm{off}_i (1-\ln n^\mathrm{off}_i)

Some special cases need to be handled separately in ``wstat``.

If :math:`n^\mathrm{on}_i = 0` but :math:`n^\mathrm{off}_i > 0` the
contribution to the log-likelihood from the energy bin :math:`i` is

.. math::
   -\ln L_i(M_s) = s_i(M_s) + n^\mathrm{off}_i \ln(\alpha_i+1).

If :math:`n^\mathrm{off}_i = 0` and
:math:`n^\mathrm{on}_i > s_i(M_s) \frac{\alpha_i + 1}{\alpha_i}`
the contribution to the log-likelihood from the energy bin :math:`i` is

.. math::
   -\ln L_i(M_s) = -\frac{s_i(M_s)}{\alpha_i} - n^\mathrm{on}_i
                   \ln\left(\frac{\alpha_i}{\alpha_i+1}\right)

However, for smaller :math:`n^\mathrm{on}_i` the value of :math:`b_i(M_s)` is
null or negative. Since a negative number of background counts is unphysical,
the number of background counts is forced to be zero. This yields the following
expression for the log-likelihood in the energy bin :math:`i`:

.. math::
   -\ln L_i(M_s) = s_i(M_s) + n^\mathrm{on}_i \left( \ln n^\mathrm{on}_i - \ln s_i(M_s) - 1 \right)

or, if also :math:`n^\mathrm{on}_i = 0`,

.. math::
   -\ln L_i(M_s) = s_i(M_s)

Forcing the number of expected background counts to zero biases the likelihood
estimator. Therefore, ``wstat`` is known to be inaccurate if there are energy
bins with zero Off counts.


