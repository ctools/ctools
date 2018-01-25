.. _howto_fermi_errors:

Determine asymmetric model parameter errors for a model component
-----------------------------------------------------------------

  .. admonition:: What you will learn

     You will learn how to use :ref:`cterror` to derive asymmetric model
     parameter errors from Fermi-LAT data using a likelihood profile method.

The :ref:`ctlike` tool computes errors for each fitted parameter value based
on the covariance matrix of the fit. This assumes a Gaussian approximation
of the uncertainties around the best fitting values. Uncertainties may however
by asymmetric.

To test the symmetry of the uncertainties you can use the :ref:`cterror`
tool. :ref:`cterror` will make use of the likelihood profile to determine
the confidence interval for each fitted parameter, allowing thus also for
asymmetric errors. You use :ref:`cterror` as follows:

.. code-block:: bash

   $ cterror
   Input event list, counts cube or observation definition XML file [events.fits] obs.xml
   Source of interest [Crab] Vela
   Input model definiton XML file [$CTOOLS/share/models/crab.xml] vela_results.xml
   Output model definiton XML file [results.xml] vela_errors.xml

The confidence intervals can then be extracted from the ``cterror.log`` log file
that was created by :ref:`cterror`. As you will see, the parameter errors
are essentially symmetric and very close to the values determined using the
covariance matrix.

.. code-block:: none

   2018-01-25T20:45:22: +================================================+
   2018-01-25T20:45:22: | Compute error for source "Vela" parameter "RA" |
   2018-01-25T20:45:22: +================================================+
   2018-01-25T20:45:22:  Confidence level ..........: 68 %
   2018-01-25T20:45:22:  Log-likelihood difference .: 0.494473
   2018-01-25T20:45:22:  Initial factor range ......: [128.819554251926, 128.851989886434]
   2018-01-25T20:45:41:  Lower parameter factor ....: 128.834
   2018-01-25T20:45:41:  Upper parameter factor ....: 128.837
   2018-01-25T20:45:41:  Error from curvature ......: 0.00162178172537938 deg
   2018-01-25T20:45:41:  Error from profile ........: 0.00169067577327553 deg
   2018-01-25T20:45:41:  Negative profile error ....: 0.00165900034895117 deg
   2018-01-25T20:45:41:  Positive profile error ....: 0.00172235119759989 deg
   2018-01-25T20:45:41:
   2018-01-25T20:45:41: +=================================================+
   2018-01-25T20:45:41: | Compute error for source "Vela" parameter "DEC" |
   2018-01-25T20:45:41: +=================================================+
   2018-01-25T20:45:41:  Confidence level ..........: 68 %
   2018-01-25T20:45:41:  Log-likelihood difference .: 0.494473
   2018-01-25T20:45:41:  Initial factor range ......: [-45.1949881923928, -45.1721741961801]
   2018-01-25T20:46:01:  Lower parameter factor ....: -45.1847
   2018-01-25T20:46:01:  Upper parameter factor ....: -45.1824
   2018-01-25T20:46:01:  Error from curvature ......: 0.00114069981063169 deg
   2018-01-25T20:46:01:  Error from profile ........: 0.00116548552428952 deg
   2018-01-25T20:46:01:  Negative profile error ....: 0.00115852324517363 deg
   2018-01-25T20:46:01:  Positive profile error ....: 0.00117244780340542 deg
   2018-01-25T20:46:01:
   2018-01-25T20:46:01: +=======================================================+
   2018-01-25T20:46:01: | Compute error for source "Vela" parameter "Prefactor" |
   2018-01-25T20:46:01: +=======================================================+
   2018-01-25T20:46:01:  Confidence level ..........: 68 %
   2018-01-25T20:46:01:  Log-likelihood difference .: 0.494473
   2018-01-25T20:46:01:  Initial factor range ......: [1e-07, 9.09784331905335]
   2018-01-25T20:46:45:  Lower parameter factor ....: 3.93674
   2018-01-25T20:46:45:  Upper parameter factor ....: 4.89558
   2018-01-25T20:46:45:  Error from curvature ......: 4.74564602380296e-10 ph/cm2/s/MeV
   2018-01-25T20:46:45:  Error from profile ........: 4.79421031311828e-10 ph/cm2/s/MeV
   2018-01-25T20:46:45:  Negative profile error ....: 4.15456324058327e-10 ph/cm2/s/MeV
   2018-01-25T20:46:45:  Positive profile error ....: 5.4338573856533e-10 ph/cm2/s/MeV
   2018-01-25T20:46:45:
   2018-01-25T20:46:45: +====================================================+
   2018-01-25T20:46:45: | Compute error for source "Vela" parameter "Index1" |
   2018-01-25T20:46:45: +====================================================+
   2018-01-25T20:46:45:  Confidence level ..........: 68 %
   2018-01-25T20:46:45:  Log-likelihood difference .: 0.494473
   2018-01-25T20:46:45:  Initial factor range ......: [1.0323933217421, 1.65461803910864]
   2018-01-25T20:47:31:  Lower parameter factor ....: 1.31084
   2018-01-25T20:47:31:  Upper parameter factor ....: 1.37275
   2018-01-25T20:47:31:  Error from curvature ......: 0.0311112358683271
   2018-01-25T20:47:31:  Error from profile ........: 0.0309517300203693
   2018-01-25T20:47:31:  Negative profile error ....: 0.032660721248488
   2018-01-25T20:47:31:  Positive profile error ....: 0.0292427387922507
   2018-01-25T20:47:31:
   2018-01-25T20:47:31: +==========================================================+
   2018-01-25T20:47:31: | Compute error for source "Vela" parameter "CutoffEnergy" |
   2018-01-25T20:47:31: +==========================================================+
   2018-01-25T20:47:31:  Confidence level ..........: 68 %
   2018-01-25T20:47:31:  Log-likelihood difference .: 0.494473
   2018-01-25T20:47:31:  Initial factor range ......: [0.001, 2.80853448949703]
   2018-01-25T20:48:14:  Lower parameter factor ....: 0.81239
   2018-01-25T20:48:14:  Upper parameter factor ....: 1.17294
   2018-01-25T20:48:14:  Error from curvature ......: 181.899981380157 MeV
   2018-01-25T20:48:14:  Error from profile ........: 180.277267463442 MeV
   2018-01-25T20:48:14:  Negative profile error ....: 177.144641591912 MeV
   2018-01-25T20:48:14:  Positive profile error ....: 183.409893334972 MeV
   2018-01-25T20:48:14:
   2018-01-25T20:48:14: +====================================================+
   2018-01-25T20:48:14: | Compute error for source "Vela" parameter "Index2" |
   2018-01-25T20:48:14: +====================================================+
   2018-01-25T20:48:14:  Confidence level ..........: 68 %
   2018-01-25T20:48:14:  Log-likelihood difference .: 0.494473
   2018-01-25T20:48:14:  Initial factor range ......: [0.308028525758507, 0.866644263529458]
   2018-01-25T20:49:01:  Lower parameter factor ....: 0.559037
   2018-01-25T20:49:01:  Upper parameter factor ....: 0.615056
   2018-01-25T20:49:01:  Error from curvature ......: 0.0279307868885475
   2018-01-25T20:49:01:  Error from profile ........: 0.0280092058458762
   2018-01-25T20:49:01:  Negative profile error ....: 0.028299015036004
   2018-01-25T20:49:01:  Positive profile error ....: 0.0277193966557485
