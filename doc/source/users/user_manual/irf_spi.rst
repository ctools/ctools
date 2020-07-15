.. _um_irf_spi:

SPI response functions
----------------------

Formulation
~~~~~~~~~~~

A SPI event is characterised by an instrument direction, spanned by the telescope
pointing direction :math:`(\alpha_{\rm spix}, \delta_{\rm spix})` and the
pseudo-detector identifier ``DETID``, and the measured energy :math:`E'`.
The SPI response is stored in the form of FITS files that can be downloaded
using

.. code-block:: bash

   $ rsync -Lzrtv isdcarc.unige.ch::arc/FTP/arc_distr/ic_tree/prod/ $REP_BASE_PROD

where ``REP_BASE_PROD`` is the target repository where all SPI data will
reside. Each SPI response FITS file contains the response functions for a
given true energy :math:`E`. The response is presented in telescope coordinates,
and gives for a given pseudo-detector identifier and true energy the effective
area as a function of sky direction :math:`(\alpha, \delta)` with respect to
the SPI telescope X axis

.. math::
   R(p',E',t'|p,E,t) = A_{\rm eff}(\alpha, \delta |{\tt DETID}, E)
