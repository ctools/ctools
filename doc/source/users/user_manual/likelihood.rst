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
   dp' \, dE' \, dt'

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
   -\ln L(M_s) = \sum_i s_i(M_s) + \alpha b_i(M_s) -
                 n^\mathrm{on}_i \ln  [s_i(M_s)+ \alpha \, b_i(M_s)] + b_i(M_s) -
                 n^\mathrm{off}_i \ln b_i(M_s) -\\
                 n^\mathrm{on}_i (1-\ln n^\mathrm{on}_i) -
                 n^\mathrm{off}_i (1-\ln n^\mathrm{off}_i)

where

.. math::
   \alpha = \frac{\int_\mathrm{on} d\Omega}{\int_\mathrm{off} d\Omega}

is the ratio between the solid angles of the On region and the Off region.
The terms in the last row are added so that :math:`-2 \ln L(M_s)` follows a
:math:`\chi^2` distribution.

Some special cases need to be handled separately in ``wstat``.
If :math:`n^\mathrm{on}_i = 0` but :math:`n^\mathrm{off}_i > 0` the
contribution to the log-likelihood from the energy bin :math:`i` is

.. math::
   -\ln L_i(M_s) = s_i(M_s) + n^\mathrm{off}_i \ln(\alpha+1).

If :math:`n^\mathrm{off}_i = 0` and
:math:`n^\mathrm{on}_i > s_i(M_s) \frac{\alpha + 1}{\alpha}`
the contribution to the log-likelihood from the energy bin :math:`i` is

.. math::
   -\ln L_i(M_s) = -\frac{s_i(M_s)}{\alpha} - n^\mathrm{on}_i
                   \ln\left(\frac{\alpha}{\alpha+1}\right)

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


Levenberg-Marquardt algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ctools uses an iterative Levenberg-Marquardt algorithm for maximum likelihood
estimation. The Levenberg-Marquardt algorithm starts with an inital guess of
the :ref:`model <um_models>` parameters :math:`a_k` and iteratively
replaces this estimate by a new estimate :math:`a_k + \Delta a_k`. The
:math:`\Delta a_k` are determined by solving

.. math::
   \sum_l \alpha_{kl} (1 + \delta_{kl} \lambda) \Delta a_l = \beta_k

where

.. math::
   \alpha_{kl} = \frac{\partial^2 (-\ln L(M))}{\partial a_k \partial a_l}

is the curvature matrix

.. math::
   \beta_k = \frac{\partial (-\ln L(M))}{\partial a_k}

is the gradient and
:math:`\delta_{kl}` is the Kronecker delta that is :math:`1` for
:math:`k=l` and :math:`0` otherwise. :math:`\lambda` is a damping parameter
that initially is set to 0.001. If a Levenberg-Marquardt iteration leads to
an increase of the log-likelihood function, :math:`\lambda` is decreased by a
factor of 10. If the log-likelihood function does not improve, :math:`\lambda`
is increased by a factor of 10 and the iteration is repeated. The iterations
stop when the log-likelihood increase is less than a small value, typically
0.005; the optimiser status is then set to ``converged``. The iterations are
also stopped if the log-likelihood function does not improve for (typically)
ten iterations; the optimiser status is then set to ``stalled``. The matrix
equation is solved using a sparse matrix Cholesky decomposition. Parameters
are constrained within their parameter limits in case they exist.

Model fitting using the Levenberg-Marquardt algorithm is implemented by
:ref:`ctlike`.


Statistical parameter errors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Statistical errors on the model parameters :math:`\delta a_k` are determined
by computing the square root of the diagonal elements of the covariance matrix
:math:`C` which is the inverse of the curvature matrix:

.. math::
   \delta a_k = \sqrt{C_{kk}}

with

.. math::
   C = [\alpha]^{-1}

Inversion of :math:`[\alpha]` is again performed using a sparse matrix Cholesky
decomposition.


Detection significance
~~~~~~~~~~~~~~~~~~~~~~

The detection significance of the source model is estimated using the so
called Test Statistic (TS) which is defined as

.. math::
   \mathrm{TS} = 2 \, (\ln L(M_b) - \ln L(M_s+M_b))

where :math:`-\ln L(M_s+M_b)` is the log-likelihood value obtained when
fitting the source and the background together to the data, and
:math:`-\ln L(M_b)` is the log-likelihood value obtained when fitting only
the background model to the data.
Under the hypothesis that the model :math:`M_b` provides a satisfactory fit
of the data, :math:`TS` follows a :math:`\chi^2_n` distribution with
:math:`n` degrees of freedom, where :math:`n` is the number of free parameters
in the source model component. Therefore

.. math::
   p = \int_\mathrm{TS}^{+\infty} \chi^2_n(x) \:\: \mathrm{d}x

gives the chance probability (p-value) that the log-likelihood improves by
:math:`TS/2` when adding the source model :math:`M_s` due to statistical
fluctuations only. For :math:`n=1` the significance in Gaussian sigma
is given by :math:`\sqrt{TS}`.


Upper flux limits
~~~~~~~~~~~~~~~~~

If gamma-ray emission from a source is not detected, an upper flux limit can
be derived by determining the flux :math:`F_\mathrm{up}` that leads to a
log-likelihood decrease of :math:`\Delta ln L` with respect to the maximum
log-likelihood estimate :math:`F_\mathrm{0}`:

.. math::
   -\ln L(F_\mathrm{up}) = -\ln L(F_\mathrm{0}) + \Delta \ln L

The log-likelihood decrease :math:`\Delta ln L` is computed from the
chance probability (p-value) using

.. math::
   \Delta \ln L = (\mathrm{erf}^{-1}(p))^2

Upper limit computation is implemented by :ref:`ctulimit`.
