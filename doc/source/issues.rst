.. _issues:

Known issues
------------

Below you will find a list of known issues.

- :ref:`Broken power law has unreliable errors <issues_bplaw>`


.. _issues_bplaw:

.. topic:: Broken power law has unreliable errors

   The broken power law spectral model has unreliable errors, specifically
   for the prefactor and the break value. Errors are in general too large,
   and this is related to the fact that the law's gradient is discontinous
   in energy. There is not very much we can do about it, it's inherent in
   the law.
