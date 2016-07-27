.. _sec_stacked:

Performing a stacked analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A stacked analysis is a binned analysis where all data from multiple
observations are stacked into a single counts cube.
The event stacking is done using the :ref:`ctbin` tool.
Instead of providing to :ref:`ctbin` an event list you should now
specify the observation definition XML file ``obs.xml`` on input.
:ref:`ctbin` will then loop over all observations and collect all events
into a single counts cube:

.. code-block:: bash

  $ ctbin
  Input event list or observation definition XML file [events.fits] obs.xml
  First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 
  Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 
  Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR] 
  Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] 
  Image scale (in degrees/pixel) [0.02] 
  Size of the X axis in pixels [200] 
  Size of the Y axis in pixels [200] 
  Algorithm for defining energy bins (FILE|LIN|LOG) [LOG] 
  Start value for first energy bin in TeV [0.1] 
  Stop value for last energy bin in TeV [100.0] 
  Number of energy bins [20] 
  Output counts cube file [cntcube.fits]

You now have a stacked counts cube ``cntcube.fits`` on disk.
Before you can use that counts cube in a maximum likelihood analysis, you have
to compute the instrument response function and the background model that is
needed for the analysis.

For the instrument response function, you have to compute the total exposure
for the stacked cube (i.e. the sum of the effective areas for each observation
multiplied by the corresponding livetimes) and an effective point spread
function (i.e. the point spread function of the different observations
weighted by the corresponding exposures).
Optionally, you can also compute an effective energy dispersion (i.e. the
energy dispersion of the different observations weighted by the corresponding
exposures).
To get these informations you use the :ref:`ctexpcube`, :ref:`ctpsfcube` and
:ref:`ctedispcube` tools:

.. code-block:: bash

  $ ctexpcube
  Input event list or observation definition XML file [NONE] obs.xml
  Calibration database [prod2] 
  Instrument response function [South_0.5h] 
  Input counts cube file to extract exposure cube definition [NONE] cntcube.fits
  Output exposure cube file [expcube.fits] 

.. code-block:: bash

  $ ctpsfcube
  Input event list or observation definition XML file [NONE] obs.xml
  Calibration database [prod2] 
  Instrument response function [South_0.5h] 
  Input counts cube file to extract PSF cube definition [NONE] 
  First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 
  Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 
  Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR] 
  Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] 
  Image scale (in degrees/pixel) [1.0] 
  Size of the X axis in pixels [10] 
  Size of the Y axis in pixels [10] 
  Lower energy limit (TeV) [0.1] 
  Upper energy limit (TeV) [100.0] 
  Number of energy bins [20] 
  Output PSF cube file [psfcube.fits] 

.. code-block:: bash

  $ ctedispcube
  Input event list or observation definition XML file [NONE] obs.xml
  Calibration database [prod2]
  Instrument response function [South_0.5h]
  Input counts cube file to extract energy dispersion cube definition [NONE]
  First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63]
  Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01]
  Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR]
  Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL]
  Image scale (in degrees/pixel) [1.0]
  Size of the X axis in pixels [10]
  Size of the Y axis in pixels [10]
  Lower energy limit (TeV) [0.1]
  Upper energy limit (TeV) [100.0]
  Number of energy bins [20]
  Output energy dispersion cube file [edispcube.fits]

.. note::

   You may have noticed that for :ref:`ctexpcube` you provided an input counts
   cube, while for the other tools you specified ``NONE``.
   By providing an input counts cube you instructed :ref:`ctexpcube` to
   extract the definition of the exposure cube from the counts cube. This is
   a convenient trick to reduce the number of user parameters that you need
   to specify. You did however not apply this trick for
   :ref:`ctpsfcube` and :ref:`ctedispcube`. In fact, the point spread function
   and energy dispersion do not vary significantly on spatial scales of 0.02°,
   and using the counts cube definition for these cubes would lead to large
   response cube files with a spatial precision that is actually not needed
   (the point spread function and energy dispersion cubes are actually
   4-dimensional data cubes, hence their size increases quickly for a large
   number of spatial pixels). Therefore, you have specified a larger image
   scale of 1° for both cubes and only a small number of 10x10 spatial pixels,
   leading to point spread function and energy dispersion cubes of modest size
   (a few MB).

You provided the ``obs.xml`` file that defines all observations on input
so that the tools know which observations were combined in the :ref:`ctbin`
run.
As final step of the analysis preparation, you need to generate a
background cube using the :ref:`ctbkgcube` tool:

.. code-block:: bash

  $ ctbkgcube
  Input event list or observation definition XML file [NONE] obs.xml
  Calibration database [prod2] 
  Instrument response function [South_0.5h] 
  Input counts cube file to extract background cube definition [NONE] cntcube.fits
  Input model XML file [NONE] $CTOOLS/share/models/crab.xml
  Output background cube file [bkgcube.fits] 
  Output model XML file [NONE] model.xml

The usage of :ref:`ctbkgcube` is very similar to that of :ref:`ctexpcube`,
yet it takes the model definition XML file as an additional input parameter.
You used here the usual ``$CTOOLS/share/models/crab.xml`` model
file that is shipped with the ctools.
:ref:`ctbkgcube` provides on output the background cube file
``bkgcube.fits`` and the model definition XML file ``model.xml`` that can
be used for further analysis.
Having a look at the ``model.xml`` file illustrates how the background
modelling works:

.. code-block:: xml

  <?xml version="1.0" encoding="UTF-8" standalone="no"?>
  <source_library title="source library">
    <source name="Crab" type="PointSource">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor" value="5.7" error="0" scale="1e-16" min="1e-07" max="1000" free="1" />
        <parameter name="Index" value="2.48" error="0" scale="-1" min="0" max="5" free="1" />
        <parameter name="PivotEnergy" value="0.3" scale="1e+06" min="0.01" max="1000" free="0" />
      </spectrum>
      <spatialModel type="SkyDirFunction">
        <parameter name="RA" value="83.6331" scale="1" min="-360" max="360" free="0" />
        <parameter name="DEC" value="22.0145" scale="1" min="-90" max="90" free="0" />
      </spatialModel>
    </source>
    <source name="BackgroundModel" type="CTACubeBackground" instrument="CTA,HESS,MAGIC,VERITAS">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor" value="1" error="0" scale="1" min="0" free="1" />
        <parameter name="Index" value="0" error="0" scale="1" min="-10" max="10" free="1" />
        <parameter name="PivotEnergy" value="1" scale="1e+06" free="0" />
      </spectrum>
    </source>
  </source_library>

The Crab source component is the same that is also present in
``$CTOOLS/share/models/crab.xml`` and is not modified.
The background component, however, has been replaced by a model of
type ``CTACubeBackground``.
This model is a 3-dimensional data cube that describes the expected 
background rate as function of spatial position and energy.
The data cube is multiplied by a power law spectrum that allows to adjust
the normalization and slope of the background spectrum in the fit.
This power law could be replaced by any spectral model that is found
as an appropriate multiplicator to the background cube.

.. note::

   There is no constraint on providing the same spatial binning or
   the same energy binning for an exposure cube, a PSF cube, an energy
   dispersion cube, a background cube and a counts cube.
   ctools interpolates internally all response cubes hence any arbitrary
   appropriate binning may be used.
   Using the same binning for the exposure cube, the background cube and
   the counts cube is only a convenience.

Now you have all files at hand to perform a stacked maximum likelihood
analysis using the :ref:`ctlike` tool:

.. code-block:: bash

  $ ctlike
  Input event list, counts cube or observation definition XML file [obs.xml] cntcube.fits
  Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits
  Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits
  Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits
  Input model XML file [$CTOOLS/share/models/crab.xml] model.xml
  Output model XML file [crab_results.xml] 

:ref:`ctlike` recognises that a counts cube should be analysed and queries
for the exposure cube, the PSF cube, and the background cube file names.
If you want to consider also the energy dispersion during the maximum likelihood
fitting you should pass the hidden ``edisp`` parameter to ctlike, and the tool
will also query of the energy dispersion cube:

.. code-block:: bash

  $ ctlike edisp=yes
  Input event list, counts cube or observation definition XML file [cntcube.fits]
  Input exposure cube file (only needed for stacked analysis) [expcube.fits]
  Input PSF cube file (only needed for stacked analysis) [psfcube.fits]
  Input background cube file (only needed for stacked analysis) [bkgcube.fits]
  Input energy dispersion cube file (only needed for stacked analysis) [NONE] edispcube.fits
  Input model XML file [model.xml]
  Output model XML file [crab_results.xml]

.. warning::

   The maximum likelihood computations including energy dispersion are
   relatively time consuming, and in many situations the impact of the
   energy dispersion on the analysis results will be very small. So make
   sure that you really need energy dispersion before you are using it.

The log file of the :ref:`ctlike` run (without energy dispersion) is shown
below.

.. code-block:: none

  2016-06-29T19:54:14: +=================================+
  2016-06-29T19:54:14: | Maximum likelihood optimisation |
  2016-06-29T19:54:14: +=================================+
  2016-06-29T19:54:15:  >Iteration   0: -logL=83633.454, Lambda=1.0e-03
  2016-06-29T19:54:15:  >Iteration   1: -logL=83561.979, Lambda=1.0e-03, delta=71.475, max(|grad|)=153.136163 [Index:7]
  2016-06-29T19:54:16:  >Iteration   2: -logL=83561.823, Lambda=1.0e-04, delta=0.156, max(|grad|)=-0.183495 [Prefactor:6]
  2016-06-29T19:54:16:  >Iteration   3: -logL=83561.823, Lambda=1.0e-05, delta=0.000, max(|grad|)=-0.003347 [Index:3]
  ...
  2016-06-29T19:54:17: +=========================================+
  2016-06-29T19:54:17: | Maximum likelihood optimisation results |
  2016-06-29T19:54:17: +=========================================+
  2016-06-29T19:54:17: === GOptimizerLM ===
  2016-06-29T19:54:17:  Optimized function value ..: 83561.823
  2016-06-29T19:54:17:  Absolute precision ........: 0.005
  2016-06-29T19:54:17:  Acceptable value decrease .: 2
  2016-06-29T19:54:17:  Optimization status .......: converged
  2016-06-29T19:54:17:  Number of parameters ......: 10
  2016-06-29T19:54:17:  Number of free parameters .: 4
  2016-06-29T19:54:17:  Number of iterations ......: 3
  2016-06-29T19:54:17:  Lambda ....................: 1e-06
  2016-06-29T19:54:17:  Maximum log likelihood ....: -83561.823
  2016-06-29T19:54:17:  Observed events  (Nobs) ...: 35946.000
  2016-06-29T19:54:17:  Predicted events (Npred) ..: 35946.000 (Nobs - Npred = 1.56502e-05)
  2016-06-29T19:54:17: === GModels ===
  2016-06-29T19:54:17:  Number of models ..........: 2
  2016-06-29T19:54:17:  Number of parameters ......: 10
  2016-06-29T19:54:17: === GModelSky ===
  2016-06-29T19:54:17:  Name ......................: Crab
  2016-06-29T19:54:17:  Instruments ...............: all
  2016-06-29T19:54:17:  Instrument scale factors ..: unity
  2016-06-29T19:54:17:  Observation identifiers ...: all
  2016-06-29T19:54:17:  Model type ................: PointSource
  2016-06-29T19:54:17:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
  2016-06-29T19:54:17:  Number of parameters ......: 6
  2016-06-29T19:54:17:  Number of spatial par's ...: 2
  2016-06-29T19:54:17:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
  2016-06-29T19:54:17:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
  2016-06-29T19:54:17:  Number of spectral par's ..: 3
  2016-06-29T19:54:17:   Prefactor ................: 5.70255e-16 +/- 7.24185e-18 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
  2016-06-29T19:54:17:   Index ....................: -2.46568 +/- 0.0110458 [-0,-5]  (free,scale=-1,gradient)
  2016-06-29T19:54:17:   PivotEnergy ..............: 300000 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2016-06-29T19:54:17:  Number of temporal par's ..: 1
  2016-06-29T19:54:17:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
  2016-06-29T19:54:17: === GCTAModelCubeBackground ===
  2016-06-29T19:54:17:  Name ......................: BackgroundModel
  2016-06-29T19:54:17:  Instruments ...............: CTA, HESS, MAGIC, VERITAS
  2016-06-29T19:54:17:  Instrument scale factors ..: unity
  2016-06-29T19:54:17:  Observation identifiers ...: all
  2016-06-29T19:54:17:  Model type ................: "PowerLaw" * "Constant"
  2016-06-29T19:54:17:  Number of parameters ......: 4
  2016-06-29T19:54:17:  Number of spectral par's ..: 3
  2016-06-29T19:54:17:   Prefactor ................: 0.964885 +/- 0.0109395 [0.01,100] ph/cm2/s/MeV (free,scale=1,gradient)
  2016-06-29T19:54:17:   Index ....................: 0.0220427 +/- 0.00686292 [-5,5]  (free,scale=1,gradient)
  2016-06-29T19:54:17:   PivotEnergy ..............: 1e+06 MeV (fixed,scale=1e+06,gradient)
  2016-06-29T19:54:17:  Number of temporal par's ..: 1
  2016-06-29T19:54:17:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
