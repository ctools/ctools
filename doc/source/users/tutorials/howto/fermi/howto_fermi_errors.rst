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

   2019-04-04T15:04:51: +================================================+
   2019-04-04T15:04:51: | Compute error for source "Vela" parameter "RA" |
   2019-04-04T15:04:51: +================================================+
   2019-04-04T15:04:51:  Confidence level ..........: 68 %
   2019-04-04T15:04:51:  Log-likelihood difference .: 0.494473
   2019-04-04T15:04:51:  Initial factor range ......: [128.819554251886, 128.851989886469]
   2019-04-04T15:05:10:  Lower parameter factor ....: 128.834
   2019-04-04T15:05:10:  Upper parameter factor ....: 128.837
   2019-04-04T15:05:10:  Error from curvature ......: 0.00162178172918532 deg
   2019-04-04T15:05:10:  Error from profile ........: 0.00169067577726878 deg
   2019-04-04T15:05:10:  Negative profile error ....: 0.00165900035287336 deg
   2019-04-04T15:05:10:  Positive profile error ....: 0.0017223512016642 deg
   2019-04-04T15:05:10:
   2019-04-04T15:05:10: +=================================================+
   2019-04-04T15:05:10: | Compute error for source "Vela" parameter "DEC" |
   2019-04-04T15:05:10: +=================================================+
   2019-04-04T15:05:10:  Confidence level ..........: 68 %
   2019-04-04T15:05:10:  Log-likelihood difference .: 0.494473
   2019-04-04T15:05:10:  Initial factor range ......: [-45.1949881925207, -45.1721741960936]
   2019-04-04T15:05:31:  Lower parameter factor ....: -45.1847
   2019-04-04T15:05:31:  Upper parameter factor ....: -45.1824
   2019-04-04T15:05:31:  Error from curvature ......: 0.00114069982135549 deg
   2019-04-04T15:05:31:  Error from profile ........: 0.0011654855352532 deg
   2019-04-04T15:05:31:  Negative profile error ....: 0.00115852325606625 deg
   2019-04-04T15:05:31:  Positive profile error ....: 0.00117244781444015 deg
   2019-04-04T15:05:31:
   2019-04-04T15:05:31: +=======================================================+
   2019-04-04T15:05:31: | Compute error for source "Vela" parameter "Prefactor" |
   2019-04-04T15:05:31: +=======================================================+
   2019-04-04T15:05:31:  Confidence level ..........: 68 %
   2019-04-04T15:05:31:  Log-likelihood difference .: 0.494473
   2019-04-04T15:05:31:  Initial factor range ......: [1e-07, 9.09784331943508]
   2019-04-04T15:06:15:  Lower parameter factor ....: 3.93674
   2019-04-04T15:06:15:  Upper parameter factor ....: 4.89558
   2019-04-04T15:06:15:  Error from curvature ......: 4.74564602405264e-10 ph/cm2/s/MeV
   2019-04-04T15:06:15:  Error from profile ........: 4.79421031332425e-10 ph/cm2/s/MeV
   2019-04-04T15:06:15:  Negative profile error ....: 4.15456324070932e-10 ph/cm2/s/MeV
   2019-04-04T15:06:15:  Positive profile error ....: 5.43385738593917e-10 ph/cm2/s/MeV
   2019-04-04T15:06:15:
   2019-04-04T15:06:15: +====================================================+
   2019-04-04T15:06:15: | Compute error for source "Vela" parameter "Index1" |
   2019-04-04T15:06:15: +====================================================+
   2019-04-04T15:06:15:  Confidence level ..........: 68 %
   2019-04-04T15:06:15:  Log-likelihood difference .: 0.494473
   2019-04-04T15:06:15:  Initial factor range ......: [1.03239332172832, 1.65461803910225]
   2019-04-04T15:07:01:  Lower parameter factor ....: 1.31084
   2019-04-04T15:07:01:  Upper parameter factor ....: 1.37275
   2019-04-04T15:07:01:  Error from curvature ......: 0.0311112358686962
   2019-04-04T15:07:01:  Error from profile ........: 0.0309517300207367
   2019-04-04T15:07:01:  Negative profile error ....: 0.0326607212488756
   2019-04-04T15:07:01:  Positive profile error ....: 0.0292427387925978
   2019-04-04T15:07:01:
   2019-04-04T15:07:01: +==========================================================+
   2019-04-04T15:07:01: | Compute error for source "Vela" parameter "CutoffEnergy" |
   2019-04-04T15:07:01: +==========================================================+
   2019-04-04T15:07:01:  Confidence level ..........: 68 %
   2019-04-04T15:07:01:  Log-likelihood difference .: 0.494473
   2019-04-04T15:07:01:  Initial factor range ......: [0.001, 2.80853448937906]
   2019-04-04T15:07:45:  Lower parameter factor ....: 0.81239
   2019-04-04T15:07:45:  Upper parameter factor ....: 1.17294
   2019-04-04T15:07:45:  Error from curvature ......: 181.899981373502 MeV
   2019-04-04T15:07:45:  Error from profile ........: 180.277267455479 MeV
   2019-04-04T15:07:45:  Negative profile error ....: 177.144641582696 MeV
   2019-04-04T15:07:45:  Positive profile error ....: 183.409893328263 MeV
   2019-04-04T15:07:45:
   2019-04-04T15:07:45: +====================================================+
   2019-04-04T15:07:45: | Compute error for source "Vela" parameter "Index2" |
   2019-04-04T15:07:45: +====================================================+
   2019-04-04T15:07:45:  Confidence level ..........: 68 %
   2019-04-04T15:07:45:  Log-likelihood difference .: 0.494473
   2019-04-04T15:07:45:  Initial factor range ......: [0.3080285257564, 0.866644263519715]
   2019-04-04T15:08:32:  Lower parameter factor ....: 0.559037
   2019-04-04T15:08:32:  Upper parameter factor ....: 0.615056
   2019-04-04T15:08:32:  Error from curvature ......: 0.0279307868881658
   2019-04-04T15:08:32:  Error from profile ........: 0.0280092058454933
   2019-04-04T15:08:32:  Negative profile error ....: 0.028299015035617
   2019-04-04T15:08:32:  Positive profile error ....: 0.0277193966553695 
