.. _howto_fermi_fitting:

Fit a model to the data
-----------------------

  .. admonition:: What you will learn

     You will learn how to use :ref:`ctlike` to fit a parametric model to
     Fermi-LAT data.

Start with creating an
:ref:`observation definition file <glossary_obsdef>`
that defines the names of the counts map file, the exposure map file, the
live time cube file and the name of the instrument response function that
should be used for the analysis:

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <observation_list title="observation library">
     <observation name="Vela" id="000001" instrument="LAT">
       <parameter name="CountsMap"    file="srcmaps.fits"/>
       <parameter name="ExposureMap"  file="expmap.fits"/>
       <parameter name="LiveTimeCube" file="ltcube.fits"/>
       <parameter name="IRF"          value="P8R2_SOURCE_V6"/>
     </observation>
   </observation_list>

This will be your input file to all tools.
Then you also need a
:ref:`model definition file <glossary_moddef>`
that describes the source and background model that you want to fit. In the
example below, a point source with a super exponentially cut-off power law is
fit at the position of the Vela pulsar on top of a Galactic diffuse and an
extragalactic isotropic source component.

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <source_library title="source library">
     <source type="PointSource" name="Vela">
       <spectrum type="SuperExponentialCutoffPowerLaw">
         <parameter name="Prefactor"    scale="1e-9"  value="1.0" min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index1"       scale="-1"    value="2.0" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="CutoffEnergy" scale="1e3"   value="1.0" min="1e-3"  max="1e+3"   free="1"/>
         <parameter name="PivotEnergy"  scale="1e3"   value="1.0" min="0.01"  max="1000.0" free="0"/>
         <parameter name="Index2"       scale="1.0"   value="0.5" min="0.1"   max="5.0"    free="1"/>
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA"  scale="1.0" value="128.84" min="-360" max="360" free="1"/>
         <parameter name="DEC" scale="1.0" value="-45.18" min="-90"  max="90"  free="1"/>
       </spatialModel>
     </source>
     <source type="DiffuseSource" name="Galactic_diffuse">
       <spectrum type="Constant">
         <parameter name="Normalization" scale="1.0" value="1.0" min="0.1" max="1000.0" free="1"/>
       </spectrum>
       <spatialModel type="DiffuseMapCube" file="gll_iem_v06.fits">
         <parameter name="Normalization" scale="1.0" value="1.0" min="0.1" max="10.0" free="0"/>
       </spatialModel>
     </source>
     <source type="DiffuseSource" name="Extragalactic_diffuse">
       <spectrum type="FileFunction" file="iso_P8R2_SOURCE_V6_v06.txt">
         <parameter name="Normalization" scale="1.0" value="1.0" min="0.0" max="1000.0" free="0"/>
       </spectrum>
       <spatialModel type="DiffuseIsotropic">
         <parameter name="Value" scale="1.0" value="1.0" min="0.0" max="10.0" free="0"/>
       </spatialModel>
     </source>
   </source_library>

.. warning::
   The names of the diffuse model components have to be identical to the
   names used when preparing the ``srcmaps.fits`` file using ``gtsrcmap``.

Now you are ready to fit the model using the :ref:`ctlike` model tool. You do
this as follows:

.. code-block:: bash

   $ ctlike
   Input event list, counts cube or observation definition XML file [events.fits] obs.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models.xml
   Output model definition XML file [crab_results.xml] vela_results.xml

:ref:`ctlike` creates on output the
:ref:`model definition file <glossary_moddef>`
``vela_results.xml`` where the parameter values were updated by their
fitted values, and where the statistical uncertainties were added using the
``error`` attribute. To investigate how the fit went you should inspect the
log file ``ctlike.log`` that was also created by :ref:`ctlike`:

.. code-block:: none

   2019-04-04T14:34:18: +=================================+
   2019-04-04T14:34:18: | Maximum likelihood optimisation |
   2019-04-04T14:34:18: +=================================+
   2019-04-04T14:34:20:  >Iteration   0: -logL=-177232.327, Lambda=1.0e-03
   2019-04-04T14:34:20:  >Iteration   1: -logL=-198931.136, Lambda=1.0e-03, delta=21698.809, step=1.0e+00, max(|grad|)=-22717.063466 [CutoffEnergy:4]
   2019-04-04T14:34:21:  >Iteration   2: -logL=-200269.888, Lambda=1.0e-04, delta=1338.753, step=1.0e+00, max(|grad|)=-67629.683813 [CutoffEnergy:4]
   2019-04-04T14:34:21:  >Iteration   3: -logL=-212351.965, Lambda=1.0e-05, delta=12082.076, step=1.0e+00, max(|grad|)=23026.802441 [Index2:6]
   2019-04-04T14:34:22:  >Iteration   4: -logL=-213954.259, Lambda=1.0e-06, delta=1602.294, step=1.0e+00, max(|grad|)=-8297.380921 [CutoffEnergy:4]
   2019-04-04T14:34:22:  >Iteration   5: -logL=-215347.973, Lambda=1.0e-07, delta=1393.714, step=1.0e+00, max(|grad|)=2746.754045 [Index2:6]
   2019-04-04T14:34:23:  >Iteration   6: -logL=-215497.807, Lambda=1.0e-08, delta=149.834, step=1.0e+00, max(|grad|)=699.848487 [Index2:6]
   2019-04-04T14:34:23:   Iteration   7: -logL=-215497.807, Lambda=1.0e-09, delta=-13.845, step=1.0e+00, max(|grad|)=-1201.982876 [Index1:3] (stalled)
   2019-04-04T14:34:24:   Iteration   8: -logL=-215497.807, Lambda=1.0e-08, delta=-13.844, step=1.0e+00, max(|grad|)=-1201.959648 [Index1:3] (stalled)
   2019-04-04T14:34:24:   Iteration   9: -logL=-215497.807, Lambda=1.0e-07, delta=-13.835, step=1.0e+00, max(|grad|)=-1201.727402 [Index1:3] (stalled)
   2019-04-04T14:34:24:   Iteration  10: -logL=-215497.807, Lambda=1.0e-06, delta=-13.744, step=1.0e+00, max(|grad|)=-1199.408718 [Index1:3] (stalled)
   2019-04-04T14:34:25:   Iteration  11: -logL=-215497.807, Lambda=1.0e-05, delta=-12.858, step=1.0e+00, max(|grad|)=-1176.594022 [Index1:3] (stalled)
   2019-04-04T14:34:26:   Iteration  12: -logL=-215497.807, Lambda=1.0e-04, delta=-6.126, step=1.0e+00, max(|grad|)=-980.833964 [Index1:3] (stalled)
   2019-04-04T14:34:26:  >Iteration  13: -logL=-215503.803, Lambda=1.0e-03, delta=5.995, step=1.0e+00, max(|grad|)=-285.964915 [Index1:3]
   2019-04-04T14:34:27:   Iteration  14: -logL=-215503.803, Lambda=1.0e-04, delta=-8.762, step=1.0e+00, max(|grad|)=-730.461327 [Index1:3] (stalled)
   2019-04-04T14:34:27:  >Iteration  15: -logL=-215504.976, Lambda=1.0e-03, delta=1.173, step=1.0e+00, max(|grad|)=-153.941337 [Index1:3]
   2019-04-04T14:34:28:   Iteration  16: -logL=-215503.593, Lambda=1.0e-04, delta=-1.383, step=1.0e+00, max(|grad|)=-327.368446 [Index1:3] (stalled)
   2019-04-04T14:34:28:  >Iteration  17: -logL=-215505.817, Lambda=1.0e-03, delta=2.225, step=1.0e+00, max(|grad|)=3.256215 [Index2:6]
   2019-04-04T14:34:29:  >Iteration  18: -logL=-215505.824, Lambda=1.0e-04, delta=0.007, step=1.0e+00, max(|grad|)=-8.190315 [CutoffEnergy:4]
   2019-04-04T14:34:29:  >Iteration  19: -logL=-215505.826, Lambda=1.0e-05, delta=0.002, step=1.0e+00, max(|grad|)=0.635264 [Index2:6]
   2019-04-04T14:34:30:
   2019-04-04T14:34:30: +=========================================+
   2019-04-04T14:34:30: | Maximum likelihood optimisation results |
   2019-04-04T14:34:30: +=========================================+
   2019-04-04T14:34:30: === GOptimizerLM ===
   2019-04-04T14:34:30:  Optimized function value ..: -215505.826
   2019-04-04T14:34:30:  Absolute precision ........: 0.005
   2019-04-04T14:34:30:  Acceptable value decrease .: 2
   2019-04-04T14:34:30:  Optimization status .......: converged
   2019-04-04T14:34:30:  Number of parameters ......: 14
   2019-04-04T14:34:30:  Number of free parameters .: 7
   2019-04-04T14:34:30:  Number of iterations ......: 19
   2019-04-04T14:34:30:  Lambda ....................: 1e-06
   2019-04-04T14:34:30:  Maximum log likelihood ....: 215505.826
   2019-04-04T14:34:30:  Observed events  (Nobs) ...: 202330.000
   2019-04-04T14:34:30:  Predicted events (Npred) ..: 202318.925 (Nobs - Npred = 11.0749567862367)
   2019-04-04T14:34:30: === GModels ===
   2019-04-04T14:34:30:  Number of models ..........: 3
   2019-04-04T14:34:30:  Number of parameters ......: 14
   2019-04-04T14:34:30: === GModelSky ===
   2019-04-04T14:34:30:  Name ......................: Vela
   2019-04-04T14:34:30:  Instruments ...............: all
   2019-04-04T14:34:30:  Observation identifiers ...: all
   2019-04-04T14:34:30:  Model type ................: PointSource
   2019-04-04T14:34:30:  Model components ..........: "PointSource" * "SuperExponentialCutoffPowerLaw" * "Constant"
   2019-04-04T14:34:30:  Number of parameters ......: 8
   2019-04-04T14:34:30:  Number of spatial par's ...: 2
   2019-04-04T14:34:30:   RA .......................: 128.835772118406 +/- 0.00162177814111823 [-360,360] deg (free,scale=1)
   2019-04-04T14:34:30:   DEC ......................: -45.1835812138759 +/- 0.00114069728372229 [-90,90] deg (free,scale=1)
   2019-04-04T14:34:30:  Number of spectral par's ..: 5
   2019-04-04T14:34:30:   Prefactor ................: 4.35032346843698e-09 +/- 4.73845090232622e-10 [1e-16,1e-06] ph/cm2/s/MeV (free,scale=1e-09,gradient)
   2019-04-04T14:34:30:   Index1 ...................: -1.34367731364905 +/- 0.0310888116010815 [-0,-5]  (free,scale=-1,gradient)
   2019-04-04T14:34:30:   CutoffEnergy .............: 990.275379021926 +/- 181.814299577731 [1,1000000] MeV (free,scale=1000,gradient)
   2019-04-04T14:34:30:   PivotEnergy ..............: 1000 [10,1000000] MeV (fixed,scale=1000,gradient)
   2019-04-04T14:34:30:   Index2 ...................: 0.58755011716699 +/- 0.027919253961318 [0.1,5]  (free,scale=1,gradient)
   2019-04-04T14:34:30:  Number of temporal par's ..: 1
   2019-04-04T14:34:30:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-04T14:34:30:  Number of scale par's .....: 0
   2019-04-04T14:34:30: === GModelSky ===
   2019-04-04T14:34:30:  Name ......................: Galactic_diffuse
   2019-04-04T14:34:30:  Instruments ...............: all
   2019-04-04T14:34:30:  Observation identifiers ...: all
   2019-04-04T14:34:30:  Model type ................: DiffuseSource
   2019-04-04T14:34:30:  Model components ..........: "DiffuseMapCube" * "Constant" * "Constant"
   2019-04-04T14:34:30:  Number of parameters ......: 3
   2019-04-04T14:34:30:  Number of spatial par's ...: 1
   2019-04-04T14:34:30:   Normalization ............: 1 [0.1,10]  (fixed,scale=1,gradient)
   2019-04-04T14:34:30:  Number of spectral par's ..: 1
   2019-04-04T14:34:30:   Normalization ............: 1.10941973066989 +/- 0.00503412002600817 [0.1,1000]  (free,scale=1,gradient)
   2019-04-04T14:34:30:  Number of temporal par's ..: 1
   2019-04-04T14:34:30:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-04T14:34:30:  Number of scale par's .....: 0
   2019-04-04T14:34:30: === GModelSky ===
   2019-04-04T14:34:30:  Name ......................: Extragalactic_diffuse
   2019-04-04T14:34:30:  Instruments ...............: all
   2019-04-04T14:34:30:  Observation identifiers ...: all
   2019-04-04T14:34:30:  Model type ................: DiffuseSource
   2019-04-04T14:34:30:  Model components ..........: "DiffuseIsotropic" * "FileFunction" * "Constant"
   2019-04-04T14:34:30:  Number of parameters ......: 3
   2019-04-04T14:34:30:  Number of spatial par's ...: 1
   2019-04-04T14:34:30:   Value ....................: 1 [0,10]  (fixed,scale=1,gradient)
   2019-04-04T14:34:30:  Number of spectral par's ..: 1
   2019-04-04T14:34:30:   Normalization ............: 1 [0,1000]  (fixed,scale=1,gradient)
   2019-04-04T14:34:30:  Number of temporal par's ..: 1
   2019-04-04T14:34:30:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-04T14:34:30:  Number of scale par's .....: 0

The fit converged after 19 iterations with spectral parameters that are
reasonably close to those found in
`Abdo et al. 2010, ApJ, 713, 154 <http://iopscience.iop.org/article/10.1088/0004-637X/713/1/154/pdf>`_
using the Fermi-LAT Science Tools.

.. note::
   The fit results can be compared to the corresponding output of gtlike
   when the tool is applied to the same test data. As you will see below,
   the fitted model parameters are basically identical.

   .. code-block:: none

      Extragalactic_diffuse:
      Normalization: 1
      Flux: 0.000149604 photons/cm^2/s

      Galactic_diffuse:
      Value: 1.11394 +/- 0.00508393
      Flux: 0.00054415 +/- 2.48325e-06 photons/cm^2/s

      Vela:
      Prefactor: 4.3344 +/- 0.65258
      Index1: 1.34634 +/- 0.0429336
      Scale: 1
      Cutoff: 1.00992 +/- 0.25512
      Index2: 0.595018 +/- 0.0398585
      TS value: 326648
      Flux: 1.0758e-05 +/- 4.33949e-08 photons/cm^2/s

      WARNING: Fit may be bad in range [100, 149.23] (MeV)
      WARNING: Fit may be bad in range [182.299, 405.972] (MeV)
      WARNING: Fit may be bad in range [6690.96, 8173.66] (MeV)
      WARNING: Fit may be bad in range [14900.5, 18202.4] (MeV)
      WARNING: Fit may be bad in range [33182.8, 40536] (MeV)
      WARNING: Fit may be bad in range [49518.7, 60491.9] (MeV)

      Total number of observed counts: 202330
      Total number of model events: 202327

      -log(Likelihood): -215606.3555

      Elapsed CPU time: 72.33

   You may also compute the Test Statistic value using :ref:`ctlike` by
   adding the ``tscalc="1"`` attribute to the
   :ref:`model definition file <glossary_moddef>`:

   .. code-block:: xml

      <?xml version="1.0" standalone="no"?>
      <source_library title="source library">
        <source type="PointSource" name="Vela" tscalc="1">
        ...
      </source_library>

   This results in

   .. code-block:: none

      2019-04-04T14:37:50:  Name ......................: Vela
      2019-04-04T14:37:50:  Instruments ...............: all
      2019-04-04T14:37:50:  Test Statistic ............: 326341.43991719

   which is also close to the value determined by gtlike.

.. note::
   Unlike gtlike, :ref:`ctlike` can also fit the source location. This has been
   done in the example above by setting the parameter attributes to ``free="1"``.

.. warning::
   ctools supports for the moment only the fitting of point sources for
   Fermi-LAT data. Other spatial shapes need to be preconvolved with the
   instrument response functions and added to the ``srcmaps.fits`` file.
