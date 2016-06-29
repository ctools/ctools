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
Before you can use that counts cube in a maximum likelihood
analysis, you have to compute the instrument response and the
background model that are needed to describe the stacked data.

For the former, you have to compute the total exposure for the stacked
cube (i.e. the sum of the effective areas for each observation multiplied
by the corresponding livetimes) and an effective point spread function
(i.e. the point spread function of the different observations weighted by
the corresponding exposures).
To get both informations you use the :ref:`ctexpcube` and 
:ref:`ctpsfcube` tools:

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

The ``obs.xml`` file has been provided on input to specify for both tools
which observations have been combined.
You further provided the counts cube on input so that :ref:`ctexpcube` can 
read the cub definition from that file and apply it to the exposure cube.
This is a trick to reduce the number of user parameters that you need to 
provide.

You do not apply this trick when using :ref:`ctpsfcube` as this 
would lead to a hugh output file owing to the fine spatial
pixelisation of the counts cube.
Such a fine binning is not needed for the PSF cube, as the PSF evolves 
only slowly over the field of view.
It is thus sufficient to compute a PSF cube with a rather coarse spatial 
binning; here you used a spatial binning of 1 degree covering a grid of
10 x 10 degrees.

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
yet it takes the model XML file as an additional input parameter.
You used here the usual ``$CTOOLS/share/models/crab.xml`` model
file that is shipped with the ctools.
:ref:`ctbkgcube` provides on output the background cube file
``bkgcube.fits`` and the model XML file ``model.xml`` that can
be used for further analysis.
Having a look at ``model.xml`` illustrates how the background
modelling works:

.. code-block:: xml

  <?xml version="1.0" encoding="UTF-8" standalone="no"?>
  <source_library title="source library">
    <source name="Crab" type="PointSource">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor" value="5.7" error="0" scale="1e-16" min="1e-07" max="1000" free="1" />
        <parameter name="Index" value="2.48" error="0" scale="-1" min="0" max="5" free="1" />
        <parameter name="Scale" value="0.3" scale="1e+06" min="0.01" max="1000" free="0" />
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
        <parameter name="Scale" value="1" scale="1e+06" free="0" />
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
   the same energy binning for an exposure cube, a PSF cube,
   a background cube and a counts cube.
   ctools interpolates internally the exposure cube, PSF cube and
   background cube values, hence any arbitrary appropriate binning may
   be used.
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
You specified here the names of the files produced by the :ref:`ctexpcube`,
the :ref:`ctpsfcube` and the :ref:`ctbkgcube` tools.
Furthermore you provided as model the ``model.xml`` file that has been
generated by the :ref:`ctbkgcube` tool.

The log file of the :ref:`ctlike` run is shown below.
Note that the spectral model that is multiplied with the background
cube has a Prefactor of 0.998 +/- 0.011 and an Index of 0.006 +/- 0.007,
indicating a very small correction to the actual spectrum of the background
cube.
Real life situations may of course require larger correction factors.

.. code-block:: xml

  2015-12-07T21:37:45: +=================================+
  2015-12-07T21:37:45: | Maximum likelihood optimisation |
  2015-12-07T21:37:45: +=================================+
  2015-12-07T21:37:46:  >Iteration   0: -logL=83748.453, Lambda=1.0e-03
  2015-12-07T21:37:47:  >Iteration   1: -logL=83736.039, Lambda=1.0e-03, delta=12.413, max(|grad|)=15.134616 [Index:3]
  2015-12-07T21:37:48:  >Iteration   2: -logL=83736.021, Lambda=1.0e-04, delta=0.018, max(|grad|)=0.108613 [Index:3]
  2015-12-07T21:37:49:  >Iteration   3: -logL=83736.021, Lambda=1.0e-05, delta=0.000, max(|grad|)=0.000570 [Index:3]
  ...
  2015-12-07T21:37:51: +=========================================+
  2015-12-07T21:37:51: | Maximum likelihood optimisation results |
  2015-12-07T21:37:51: +=========================================+
  2015-12-07T21:37:51: === GOptimizerLM ===
  2015-12-07T21:37:51:  Optimized function value ..: 83736.021
  2015-12-07T21:37:51:  Absolute precision ........: 0.005
  2015-12-07T21:37:51:  Acceptable value decrease .: 2
  2015-12-07T21:37:51:  Optimization status .......: converged
  2015-12-07T21:37:51:  Number of parameters ......: 10
  2015-12-07T21:37:51:  Number of free parameters .: 4
  2015-12-07T21:37:51:  Number of iterations ......: 3
  2015-12-07T21:37:51:  Lambda ....................: 1e-06
  2015-12-07T21:37:51:  Maximum log likelihood ....: -83736.021
  2015-12-07T21:37:51:  Observed events  (Nobs) ...: 35198.000
  2015-12-07T21:37:51:  Predicted events (Npred) ..: 35198.000 (Nobs - Npred = 7.6773e-07)
  2015-12-07T21:37:51: === GModels ===
  2015-12-07T21:37:51:  Number of models ..........: 2
  2015-12-07T21:37:51:  Number of parameters ......: 10
  2015-12-07T21:37:51: === GModelSky ===
  2015-12-07T21:37:51:  Name ......................: Crab
  2015-12-07T21:37:51:  Instruments ...............: all
  2015-12-07T21:37:51:  Instrument scale factors ..: unity
  2015-12-07T21:37:51:  Observation identifiers ...: all
  2015-12-07T21:37:51:  Model type ................: PointSource
  2015-12-07T21:37:51:  Model components ..........: "SkyDirFunction" * "PowerLaw" * "Constant"
  2015-12-07T21:37:51:  Number of parameters ......: 6
  2015-12-07T21:37:51:  Number of spatial par's ...: 2
  2015-12-07T21:37:51:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
  2015-12-07T21:37:51:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
  2015-12-07T21:37:51:  Number of spectral par's ..: 3
  2015-12-07T21:37:51:   Prefactor ................: 5.75289e-16 +/- 7.24749e-18 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
  2015-12-07T21:37:51:   Index ....................: -2.53122 +/- 0.0113068 [-0,-5]  (free,scale=-1,gradient)
  2015-12-07T21:37:51:   PivotEnergy ..............: 300000 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2015-12-07T21:37:51:  Number of temporal par's ..: 1
  2015-12-07T21:37:51:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
  2015-12-07T21:37:51: === GCTAModelCubeBackground ===
  2015-12-07T21:37:51:  Name ......................: BackgroundModel
  2015-12-07T21:37:51:  Instruments ...............: CTA, HESS, MAGIC, VERITAS
  2015-12-07T21:37:51:  Instrument scale factors ..: unity
  2015-12-07T21:37:51:  Observation identifiers ...: all
  2015-12-07T21:37:51:  Model type ................: "PowerLaw" * "Constant"
  2015-12-07T21:37:51:  Number of parameters ......: 4
  2015-12-07T21:37:51:  Number of spectral par's ..: 3
  2015-12-07T21:37:51:   Prefactor ................: 0.998055 +/- 0.0114979 [0.01,100] ph/cm2/s/MeV (free,scale=1,gradient)
  2015-12-07T21:37:51:   Index ....................: 0.00648796 +/- 0.00697365 [-5,5]  (free,scale=1,gradient)
  2015-12-07T21:37:51:   PivotEnergy ..............: 1e+06 MeV (fixed,scale=1e+06,gradient)
  2015-12-07T21:37:51:  Number of temporal par's ..: 1
  2015-12-07T21:37:51:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
