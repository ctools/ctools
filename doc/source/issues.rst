.. _issues:

Known issues
------------

Below you will find a list of known issues.

- :ref:`Broken power law has unreliable errors <issue_bplaw>`
- :ref:`Errors become unreliable when fitting the pivot energy <issue_pivot>`
- :ref:`Binned analysis is biased when using coarse binning <issue_binned>`
- :ref:`Shell model is biased when width is comparable to angular resolution <issue_shell>`


.. _issue_bplaw:

.. topic:: Broken power law has unreliable errors

   The broken power law spectral model has unreliable errors, specifically
   for the prefactor and the break value. Errors are in general too large,
   and this is related to the fact that the law's gradient is discontinuous
   in energy. There is not very much we can do about it, it's inherent in
   the law.

.. _issue_pivot:

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

.. _issue_binned:

.. topic:: Binned analysis is biased when using coarse binning

   When performing a binned or stacked analysis you should make sure
   that the spatial and spectral binning is sufficiently fine grained.
   The spatial binning should be better than the best angular resolution
   over the energy range of interest. Use a typical value of 0.02 degrees
   per pixel for the spatial binning and 10 bins per decade for the
   spectral binning. If the binning is too coarse, the spectral parameters 
   that are fitted will be biased.

.. _issue_shell:

.. topic:: Shell model is biased when width is comparable to angular resolution

   When the width of the shell model becomes comparable to or smaller
   than the angular resolution, the shell width tends to be overestimated
   while the shell radius tends to be underestimated.
   The fitted shell width and radius should thus not be overinterpreted
   when the width is close to the angular resolution of CTA.

