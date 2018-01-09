.. _sec_tips:

Helpful tips
============

This page summarises things that you **definitely should consider** when using
the ctools software. Make sure that you read at least this section.

- :ref:`Always inspect your fit residuals <tip_residuals>`
- :ref:`Make sure your binning is sufficiently fine grained <tip_binned>`
- :ref:`Only use energy dispersion when you really need it <tip_edisp>`
- :ref:`Only compute the Test Statistics when you really need it <tip_ts>`
- :ref:`Fit extended source with a radial disk model <tip_disk>`
- :ref:`Fix the spectral parameters of a source to compute an upper limit <tip_ulimit>`



.. _tip_residuals:

.. topic:: Always inspect your fit residuals

   **Never** trust the values of a :ref:`ctlike` model fit without having
   inspected the fit residuals. Residuals in the region of your source of
   interest should be flat.


.. _tip_binned:

.. topic:: Make sure your binning is sufficiently fine grained

   If you use binned analysis, make sure that the bin size is sufficiently fine
   grained, since the model is sampled at the bin centres. A spatial binning of
   0.02 degrees should be okay. Use at least 10 bins per decade for the energy
   binning, you may consider using more bins below 100 GeV when you analyse
   data near the energy threshold.


.. _tip_edisp:

.. topic:: Only use energy dispersion when you really need it

   You should only use the energy dispersion when you are sure that it
   impacts significantly the analysis results, or when you want to produce
   final results for a publication. Taking into account energy dispersion
   adds a further dimension to the analysis which considerably slows down the
   computations. In many cases the impact of the energy dispersion is small
   and can in a first order be neglected.


.. _tip_ts:

.. topic:: Only compute the Test Statistics when you really need it

   The Test Statistic value is useful for computing the detection significance
   of a source, but for many sources in your region of interest you may actually
   not be interested in that value. Computation of the Test Statistic value
   takes some extra time, so do not request the computation if you don't really
   need it.


.. _tip_disk:

.. topic:: Fit extended sources with a radial disk model

   The fit of extended sources is best done with a radial disk model since the
   computations for this model are very fast (model type: ``RadialDisk``). Only
   if the radial disk model is a bad representation of the spatial distribution
   or to compute final values for a publication you should switch to other
   extended models, such as the ``RadialGaussian`` or ``RadialShell`` models.
   Note that the computations for elliptical models are even slower,
   hence use these models only when you really need them, and prefer the
   elliptical disk over the elliptical Gaussian model.


.. _tip_ulimit:

.. topic:: Fix the spectral parameters of a source to compute an upper limit

   When computing an upper limit using :ref:`ctulimit`, make sure that the
   spectral parameters such as index, cutoff energy, etc. are all fixed. Only
   the prefactor or integrated flux should be left free. Computation of an
   upper limit requires to make a hypothesis on the underlying spectral energy
   density distribution, hence you have to fix the parameters of this
   hypothesis. Otherwise your upper limits are probably rather meaningless.
