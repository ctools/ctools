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

.. code-block:: bash

   2017-08-26T14:16:00: +=======================================================+
   2017-08-26T14:16:00: | Compute error for source "Vela" parameter "Prefactor" |
   2017-08-26T14:16:00: +=======================================================+
   2017-08-26T14:16:00:  Confidence level ..........: 68 %
   2017-08-26T14:16:00:  Log-likelihood difference .: 0.494473
   2017-08-26T14:16:00:  Initial factor range ......: [1e-07, 9.10916268820218]
   2017-08-26T14:16:20:  Lower parameter factor ....: 3.93549
   2017-08-26T14:16:20:  Upper parameter factor ....: 4.89285
   2017-08-26T14:16:20:  Error from curvature ......: 4.75365104151534e-10 ph/cm2/s/MeV
   2017-08-26T14:16:20:  Error from profile ........: 4.78682555317737e-10 ph/cm2/s/MeV
   2017-08-26T14:16:20:  Negative profile error ....: 4.20026137925122e-10 ph/cm2/s/MeV
   2017-08-26T14:16:20:  Positive profile error ....: 5.37338972710351e-10 ph/cm2/s/MeV
   2017-08-26T14:16:20:
   2017-08-26T14:16:20: +====================================================+
   2017-08-26T14:16:20: | Compute error for source "Vela" parameter "Index1" |
   2017-08-26T14:16:20: +====================================================+
   2017-08-26T14:16:20:  Confidence level ..........: 68 %
   2017-08-26T14:16:20:  Log-likelihood difference .: 0.494473
   2017-08-26T14:16:20:  Initial factor range ......: [1.03221907174349, 1.65479046234989]
   2017-08-26T14:16:42:  Lower parameter factor ....: 1.31121
   2017-08-26T14:16:42:  Upper parameter factor ....: 1.37307
   2017-08-26T14:16:42:  Error from curvature ......: 0.03112856953032
   2017-08-26T14:16:42:  Error from profile ........: 0.0309309760713872
   2017-08-26T14:16:42:  Negative profile error ....: 0.0322989307870749
   2017-08-26T14:16:42:  Positive profile error ....: 0.0295630213556994
   2017-08-26T14:16:42:
   2017-08-26T14:16:42: +==========================================================+
   2017-08-26T14:16:42: | Compute error for source "Vela" parameter "CutoffEnergy" |
   2017-08-26T14:16:42: +==========================================================+
   2017-08-26T14:16:42:  Confidence level ..........: 68 %
   2017-08-26T14:16:42:  Log-likelihood difference .: 0.494473
   2017-08-26T14:16:42:  Initial factor range ......: [0.001, 2.8071605417966]
   2017-08-26T14:17:02:  Lower parameter factor ....: 0.813181
   2017-08-26T14:17:02:  Upper parameter factor ....: 1.17404
   2017-08-26T14:17:02:  Error from curvature ......: 181.871842544747 MeV
   2017-08-26T14:17:02:  Error from profile ........: 180.431484426618 MeV
   2017-08-26T14:17:02:  Negative profile error ....: 175.261332662553 MeV
   2017-08-26T14:17:02:  Positive profile error ....: 185.601636190684 MeV
   2017-08-26T14:17:02:
   2017-08-26T14:17:02: +====================================================+
   2017-08-26T14:17:02: | Compute error for source "Vela" parameter "Index2" |
   2017-08-26T14:17:02: +====================================================+
   2017-08-26T14:17:02:  Confidence level ..........: 68 %
   2017-08-26T14:17:02:  Log-likelihood difference .: 0.494473
   2017-08-26T14:17:02:  Initial factor range ......: [0.307804018748382, 0.866425033775455]
   2017-08-26T14:17:24:  Lower parameter factor ....: 0.559156
   2017-08-26T14:17:24:  Upper parameter factor ....: 0.615141
   2017-08-26T14:17:24:  Error from curvature ......: 0.0279310507513537
   2017-08-26T14:17:24:  Error from profile ........: 0.0279924226890397
   2017-08-26T14:17:24:  Negative profile error ....: 0.027958327168103
   2017-08-26T14:17:24:  Positive profile error ....: 0.0280265182099764
