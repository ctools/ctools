.. _start_edisp:

Taking the energy dispersion into account
-----------------------------------------

  .. admonition:: What you will learn

     You will learn how to **take the energy dispersion (or energy redistribution)
     into account**.

     Although the effect of the energy dispersion is often neglegible there
     may be cases where you want to consider energy dispersion in an analysis,
     for example if you are analysing the data down to very low energies.

  .. warning::
     Energy dispersion is fully implemented in ctools but is computationally
     intensive. So be aware that the **tools and scripts will take a substantial
     amount of computing time if energy dispersion is considered**.

You may not have recognised, but the examples you have exercised so far
have neglected the impact of the energy dispersion on the analysis. In reality,
however, the reconstructed event energy will differ from the true photon energy,
and this effect will become particularily important at low energies. There are
therefore cases where you want to take the energy dispersion into account.

To simulate events taking the energy dispersion into account you run the
:ref:`ctobssim` tool with ``edisp=yes`` parameter:

.. code-block:: bash

   $ ctobssim edisp=yes
   RA of pointing (degrees) (0-360) [83.63]
   Dec of pointing (degrees) (-90-90) [22.01]
   Radius of FOV (degrees) (0-180) [5.0]
   Start time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:00:00]
   Stop time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:30:00] 2020-01-01T01:00:00
   Lower energy limit (TeV) [0.1] 0.030
   Upper energy limit (TeV) [100.0] 200.0
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input model definition XML file [$CTOOLS/share/models/crab.xml]
   Output event data file or observation definition XML file [events.fits] events_edisp.fits

You then select events from the simulated data as before using the
:ref:`ctselect` tool:

.. code-block:: bash

   $ ctselect
   Input event list or observation definition XML file [events.fits] events_edisp.fits
   RA for ROI centre (degrees) (0-360) [83.63]
   Dec for ROI centre (degrees) (-90-90) [22.01]
   Radius of ROI (degrees) (0-180) [3.0]
   Start time (UTC string, JD, MJD or MET in seconds) [NONE] 2020-01-01T00:10:00
   Stop time (UTC string, JD, MJD or MET in seconds) [NONE] 2020-01-01T00:40:00
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Output event list or observation definition XML file [selected_events.fits] selected_events_edisp.fits

Then you bin as before the selected events into a counts cube using the
:ref:`ctbin` tool:

.. code-block:: bash

   $ ctbin
   Input event list or observation definition XML file [selected_events.fits] selected_events_edisp.fits
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63]
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01]
   Projection method (AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL]
   Image scale (in degrees/pixel) [0.02]
   Size of the X axis in pixels [200]
   Size of the Y axis in pixels [200]
   Algorithm for defining energy bins (FILE|LIN|LOG) [LOG]
   Start value for first energy bin in TeV [0.1]
   Stop value for last energy bin in TeV [100.0]
   Number of energy bins (1-200) [20]
   Output counts cube file [cntcube.fits] cntcube_edisp.fits

As next step you need to compute the
:ref:`energy dispersion cube <glossary_edispcube>`
using the :ref:`ctexpcube` tool. You run the tool as follows:

.. code-block:: bash

   $ ctedispcube
   Input event list or observation definition XML file [NONE] selected_events_edisp.fits
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input counts cube file to extract energy dispersion cube definition [NONE]
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63]
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01]
   Projection method (AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL]
   Image scale (in degrees/pixel) [1.0]
   Size of the X axis in pixels [10]
   Size of the Y axis in pixels [10]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins [20]
   Output energy dispersion cube file [edispcube.fits]

Now you are ready to perform a binned maximum likelihood analysis taking the
energy dispersion into account. You do this by running the :ref:`ctlike` tool
with the ``edisp=yes`` parameter. The :ref:`ctlike` tool will now query for the
energy dispersion cube:

.. code-block:: bash

   $ ctlike edisp=yes
   Input event list, counts cube or observation definition XML file [selected_events.fits] cntcube_edisp.fits
   Input exposure cube file [expcube.fits]
   Input PSF cube file [psfcube.fits]
   Input energy dispersion cube file [NONE] edispcube.fits
   Input background cube file [bkgcube.fits]
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models.xml
   Output model definition XML file [crab_results.xml] crab_results_edisp.xml

You can also perform an unbinned maximum likelihood analysis taking the energy
dispersion into account. In that case the energy dispersion information will be
directly determined from the
:ref:`instrument response functions <glossary_irf>`
and no energy dispersion cube is required:

.. code-block:: bash

   $ ctlike edisp=yes
   Input event list, counts cube or observation definition XML file [cntcube_edisp.fits] selected_events_edisp.fits
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input model definition XML file [models.xml] $CTOOLS/share/models/crab.xml
   Output model definition XML file [crab_results_edisp.xml] 
