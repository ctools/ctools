.. _howto_xspec:

How to perform spectral fitting using Xspec?
--------------------------------------------

  .. admonition:: What you will learn

     You will learn how to **use Xspec to perform a spectral fitting of On/Off
     data**.

  .. warning:: You can reproduce this tutorial only if you have
     `Xspec <https://heasarc.nasa.gov/xanadu/xspec/>`_ installed.
     You may
     `download Xspec here <https://heasarc.gsfc.nasa.gov/lheasoft/download.html>`_.

If you are familiar with X-ray data analysis, or if you want to use a specific
spectral model that is available in the
`Xspec <https://heasarc.nasa.gov/xanadu/xspec/>`_
package, you can use that package to fit CTA data that were prepared for an
On/Off analysis.

To illustrate how you can do this, let's simulate an observation of the Crab
nebula with a pointing direction that is offset by 1 deg. We also use energy
dispersion in the simulation as an On/Off analysis always

.. code-block:: bash

   $ ctobssim edisp=yes
   RA of pointing (degrees) (0-360) [83.63]
   Dec of pointing (degrees) (-90-90) [22.51] 21.01
   Radius of FOV (degrees) (0-180) [5.0]
   Start time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:00:00]
   Stop time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:30:00]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input model definition XML file [$CTOOLS/share/models/crab.xml]
   Output event data file or observation definition XML file [events.fits]

No use the :ref:`csphagen` script to generate the On/Off data using a On
region radius of 0.2 deg and the reflected region model to determine the
Off counts spectrum:

.. code-block:: bash

   $ csphagen
   Input event list or observation definition XML file [obs.xml] events.fits
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input model definition XML file (if NONE, use point source) [NONE]
   Algorithm for defining energy bins (FILE|LIN|LOG) [LOG]
   Start value for first energy bin in TeV [0.1]
   Stop value for last energy bin in TeV [100.0]
   Number of energy bins [120] 20
   Stack multiple observations into single PHA, ARF and RMF files? [no]
   Output observation definition XML file [onoff_obs.xml]
   Output model definition XML file [onoff_model.xml]
   Method for background estimation (REFLECTED|CUSTOM) [REFLECTED]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL]
   Right Ascension of source region centre (deg) (0-360) [83.63]
   Declination of source region centre (deg) (-90-90) [22.01]
   Radius of source region circle (deg) (0-180) [0.2]

Now you can start
`Xspec <https://heasarc.nasa.gov/xanadu/xspec/>`_
for model fitting.

.. code-block:: bash

   $ xspec

                  XSPEC version: 12.10.1
        Build Date/Time: Mon Jan  7 12:38:07 2019

   XSPEC12>

First load the On file ``onoff_pha_on.fits`` using the ``data`` command:

.. code-block:: bash

   XSPEC12>data onoff_pha_on.fits
   ***Warning: Data file onoff_pha_on.fits has both POISSERR key set to 'true' and a STAT_ERR column.
      XSPEC will assume Poisson errors.
   ***Warning: Background file onoff_pha_off.fits has both POISSERR key set to 'true' and a STAT_ERR column.
      XSPEC will use assume Poisson errors.

   1 spectrum  in use
 
   Spectral Data File: onoff_pha_on.fits  Spectrum 1
   Net count rate (cts/s) for Spectrum:1  1.619e+00 +/- 3.135e-02 (93.9 % total)
    Assigned to Data Group 1 and Plot Group 1
     Noticed Channels:  1-20
     Telescope: CTA Instrument: PROD2  Channel Type: PI
     Exposure Time: 1764 sec
    Using fit statistic: chi
    Using test statistic: chi
    Using Background File                onoff_pha_off.fits
     Background Exposure Time: 1764 sec
    Using Response (RMF) File            onoff_rmf.fits for Source 1
    Using Auxiliary Response (ARF) File  onoff_arf.fits

This will load the On ``PHA`` file as well as the Off ``PHA``file, the ``ARF``
file and the ``RMF`` file. While you specified the On file on the command line,
the file names of the other files were extracted from ``backfile``, ``ancrfile``
and ``respfile`` header keywords in the On file that were setup by the
:ref:`csphagen` script.

Now set the fit statistic to ``cstat`` using

.. code-block:: bash

   XSPEC12>statistic cstat
   Default fit statistic is set to: C-Statistic
      This will apply to all current and newly loaded spectra.

and set the spectral model as well as the initial parameters of the spectral
model using

.. code-block:: bash

   XSPEC12>model powerlaw

   Input parameter value, delta, min, bot, top, and max values for ...
                 1       0.01(      0.01)         -3         -2          9         10
   1:powerlaw:PhoIndex>2.5
                 1       0.01(      0.01)          0          0      1e+20      1e+24
   2:powerlaw:norm>500.0

   ========================================================================
   Model powerlaw<1> Source No.: 1   Active/On
   Model Model Component  Parameter  Unit     Value
    par  comp
      1    1   powerlaw   PhoIndex            2.50000      +/-  0.0
      2    1   powerlaw   norm                500.000      +/-  0.0
   ________________________________________________________________________


   Fit statistic : C-Statistic =         727.46 using 20 PHA bins and 18 degrees of freedom.

   Test statistic : Chi-Squared =         548.66 using 20 PHA bins.
    Reduced chi-squared =         30.481 for     18 degrees of freedom
    Null hypothesis probability =  5.939476e-105

   ***Warning: Chi-square may not be valid due to bins with zero variance
               in spectrum number(s): 1

    Current data and model not fit yet.

Note that
`Xspec <https://heasarc.nasa.gov/xanadu/xspec/>`_
energies are in keV, and the power law normalization is taken at a reference
energy of 1 keV.

Finally, you are ready to do the spectral fitting using the ``fit``
command:

.. code-block:: bash

   XSPEC12>fit
                                   Parameters
   C-Statistic  |beta|/N    Lvl    1:PhoIndex        2:norm
   96.8692      11494.2      -3       2.49807       711.693
   12.0065      4626.2       -4       2.49140       735.585
   10.9389      553.368      -5       2.48822       701.803
   10.8978      108.533      -6       2.48774       697.592
   10.8974      11.7189      -7       2.48768       697.003
   ==============================
    Variances and Principal Axes
                    1        2
    8.7349E-07| -1.0000   0.0001
    4.6499E+04| -0.0001  -1.0000
   ------------------------------

   ========================
     Covariance Matrix
           1           2
      2.228e-04   3.212e+00
      3.212e+00   4.650e+04
   ------------------------

   ========================================================================
   Model powerlaw<1> Source No.: 1   Active/On
   Model Model Component  Parameter  Unit     Value
    par  comp
      1    1   powerlaw   PhoIndex            2.48768      +/-  1.49258E-02
      2    1   powerlaw   norm                697.003      +/-  215.637
   ________________________________________________________________________


   Fit statistic : C-Statistic =          10.90 using 20 PHA bins and 18 degrees of freedom.

   Test statistic : Chi-Squared =         126.30 using 20 PHA bins.
    Reduced chi-squared =         7.0169 for     18 degrees of freedom
    Null hypothesis probability =   2.683238e-18

   ***Warning: Chi-square may not be valid due to bins with zero variance
               in spectrum number(s): 1

   XSPEC12>

The fit results are compared in the table below to the simulated true values
and the values obtained using :ref:`ctlike` using the same spectral model.
The
`Xspec <https://heasarc.nasa.gov/xanadu/xspec/>`_
results are compatible with the simulated values and very close to the values
obtained using :ref:`ctlike`.

 +-----------+-------+-----------------+-----------------+
 | Parameter | Truth |      Xspec      |     ctlike      |
 +===========+=======+=================+=================+
 | Prefactor | 601.4 | 697.0 +/- 215.6 | 697.5 +/- 234.5 |
 +-----------+-------+-----------------+-----------------+
 | Index     | 2.48  | 2.488 +/- 0.015 | 2.488 +/- 0.017 |
 +-----------+-------+-----------------+-----------------+
