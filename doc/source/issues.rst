.. _issues:

Known issues
------------

Below you will find a list of known issues.

- :ref:`Broken power law has unreliable errors <issues_bplaw>`
- :ref:`Errors become unreliable when fitting the pivot energy <issues_pivot>`
- :ref:`Binned analysis is biased when using coarse binning <issues_binned>`


.. _issues_bplaw:

.. topic:: Broken power law has unreliable errors

   The broken power law spectral model has unreliable errors, specifically
   for the prefactor and the break value. Errors are in general too large,
   and this is related to the fact that the law's gradient is discontinous
   in energy. There is not very much we can do about it, it's inherent in
   the law.

.. _issues_pivot:

.. topic:: Errors become unreliable when fitting the pivot energy

   The spectral ``PowerLaw``, ``ExpCutoff`` and ``LogParabola`` models
   have a pivot energy, specified by the ``Scale`` parameter, and this
   pivot energy can not be determined in a fit together with the other
   model parameters. The reason is that the pivot energy is not an
   independent parameter of these models, and hence when all other
   spectral parameters are free, the pivot energy is unconstrained.
   So please make sure that the pivot energy is fixed, or fix other
   parameters of the model to assure non-degeneracy of the free
   parameters.

.. _issues_binned:

.. topic:: Binned analysis is biased when using coarse binning

   When performing a binned or stacked analysis you should make sure
   that the spatial and spectral binning is sufficiently fine grained.
   The spatial binning should be better than the best angular resolution
   over the energy range of interest. Use a typical value of 0.02 degrees.
   For the spectral binning, use at least 5 bins per decade.
   If the binning is too coarse, the spectral parameters will be biased.
   For example, a spectral hardening is observed when using a coarser
   spatial binning.
