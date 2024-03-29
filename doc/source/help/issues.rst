.. _issues:

Known issues
------------

Below you will find a list of known ctools issues. You may also check the
list of
`known GammaLib issues <http://cta.irap.omp.eu/gammalib/doc/html/issues.html>`_.

**Installation issues**

- :ref:`Certificate problem when cloning from Git <issue_ssl>`
- :ref:`ctools does not compile against conda Python <conda_python>`
- :ref:`Python module does not work <issue_python>`
- :ref:`Installing on Solaris <issue_solaris>`
- :ref:`Installing on OpenSolaris <issue_opensolaris>`
- :ref:`Conda installation on Ubuntu 18.04 <issue_ubuntu1804>`

**Analysis issues**

- :ref:`Binned analysis is biased for coarse binning <issue_binned>`

**Model fitting issues**

- :ref:`Pivot energy should be comprised in energy range of fitted data <issue_interval>`
- :ref:`Errors become unreliable when fitting the pivot energy <issue_pivot>`
- :ref:`Broken power law has unreliable errors <issue_bplaw>`
- :ref:`Shell model is biased when width is comparable to angular resolution <issue_shell>`
- :ref:`Elliptical Gaussian model fit converges slowly <issue_egauss>`

.. _installation_issues:

Installation issues
^^^^^^^^^^^^^^^^^^^

.. _issue_ssl:

.. topic:: Certificate problem when cloning from Git

   When cloning ctools from git you may encounter an SSL certificate 
   problem. This is related to the usage of the https protocol and can
   be resolved by setting the ``GIT_SSL_NO_VERIFY`` environment variable
   to true:

   .. code-block:: bash

      export GIT_SSL_NO_VERIFY=true

.. _conda_python:

.. topic:: ctools does not compile against conda Python

   Trying to compile ctools against conda Python may fail due to
   incompatibility issues. If you'd like to compile ctools against conda
   Python, make sure that gcc, swig, and cfitsio are all installed via
   anaconda

   .. code-block:: bash

      $ conda install gcc swig
      $ conda install -c conda-forge cfitsio

.. _issue_python:

.. topic:: Python module does not work

   ctools include a Python module that is built from so called wrapper 
   files that are autogenerated using the `swig <http://www.swig.org/>`_
   tool. These wrapper files are shipped with a ctools release, but if
   you use the code from git you need `swig <http://www.swig.org/>`_
   to generate the wrapper files during the build step. In any case,
   to compile the Python module ctools need the ``Python.h`` header file
   which may not necessarily be installed on your system. Check the output
   of ``./configure`` to examine the configuration that ctools has
   detected. You may see the following::

   * Python                       (yes)
   * Python.h                     (yes)
   * swig                         (yes)
   * Python wrappers              (yes)

   Recall, if the wrappers exist you do not need `swig <http://www.swig.org/>`_,
   but if the wrappers don't exist you need `swig <http://www.swig.org/>`_.
   If the ``Python.h`` header file does not exist then install the Python
   development package.

.. _issue_solaris:

.. topic:: Installing on Solaris

   Although ctools build on Solaris using the Sun compiler, there are
   problems with global symbols in shared libraries and exception catching,
   which prevents the FITS interface to work correctly. ctools have however
   been built and tested successfully using the GNU compiler, and this is
   the only build method that is currently supported. Problems have also
   been encountered when compiling cfitsio versions more recent than 3.250.
   The problems have been reported to the cfitsio developer team, and are
   likely to be solved in the future. For the time being, it is recommended
   to use cfitsio version 3.250 on Solaris.

.. _issue_opensolaris:

.. topic:: Installing on OpenSolaris

   On OpenSolaris, the same problems concerning the SunStudio compiler
   occur as for Solaris, and also here, the GNU compiler is the recommended
   tool to build ctools. Also here, cfitsio version 3.250 is the recommended
   library as more recent version feature relocation problems. ctools have
   been tested using gcc 4.3.2 on OpenSolaris 2009.06. Make sure to create
   the symbolic links

   .. code-block:: csh

      $ ln -s /usr/bin/gcc4.3.2 /usr/bin/gcc
      $ ln -s /usr/bin/g++4.3.2 /usr/bin/g++

   which are not there by default. This avoids warnings during compilation.

.. _issue_ubuntu1804:

.. topic:: Conda installation on Ubuntu 18.04

   According to `this thread <https://github.com/ContinuumIO/anaconda-issues/issues/11371>`_
   conda does not successfully work on Ubuntu 18.04 and consequently you will
   also encounter problems when installing ctools on Ubuntu 18.04 via conda.
   If you want to use conda on Ubuntu you need to upgrade to a newer version.
   Successful operations were reported, for example, for Ubuntu 20.04.


Analysis issues
^^^^^^^^^^^^^^^

.. _issue_binned:

.. topic:: Binned analysis is biased for coarse binning

   When performing a binned or stacked analysis you should make sure
   that the spatial and spectral binning is sufficiently fine grained.
   The spatial binning should be better than the best angular resolution
   over the energy range of interest. Use a typical value of 0.02 degrees
   per pixel for the spatial binning and 10 bins per decade for the
   spectral binning. If the binning is too coarse, the spectral parameters 
   that are fitted will be biased.


Model fitting issues
^^^^^^^^^^^^^^^^^^^^

.. _issue_interval:

.. topic:: Pivot energy should be comprised in energy range of fitted data

   The pivot energy of a spectral model, such as for example a power law model,
   should be comprised within the energy range of the fitted data,
   otherwise some fit instabilities may occur.

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

.. _issue_bplaw:

.. topic:: Broken power law has unreliable errors

   The broken power law spectral model has unreliable errors, specifically
   for the prefactor and the break value. Errors are in general too large,
   and this is related to the fact that the law's gradient is discontinuous
   in energy. There is not very much we can do about it, it's inherent in
   the law.

.. _issue_shell:

.. topic:: Shell model is biased when width is comparable to angular resolution

   When the width of the shell model becomes comparable to or smaller
   than the angular resolution, the shell width tends to be overestimated
   while the shell radius tends to be underestimated.
   The fitted shell width and radius should thus not be overinterpreted
   when the width is close to the angular resolution of CTA.

.. _issue_egauss:

.. topic:: Elliptical Gaussian model fit converges slowly

   The convergence of the elliptical Gaussian model can be slow and
   in some situations requires of the order of 20 iterations before
   the fit terminates. Nevertheless, the numerical accuracy of the model
   fitting results are satisfactory.
