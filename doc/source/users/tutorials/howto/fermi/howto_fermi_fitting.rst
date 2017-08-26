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
         <parameter name="RA"  scale="1.0" value="128.84" min="-360" max="360" free="0"/>
         <parameter name="DEC" scale="1.0" value="-45.18" min="-90"  max="90"  free="0"/>
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
..

  .. warning::

     The names of the diffuse model components have to be identical to the
     names used when preparing the ``srcmaps.fits`` file using ``gtsrcmap``.

     **Positions of point sources can so far not be fitted for a Fermi-LAT
     analysis**, hence the corresponding parameters are fixed.

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

.. code-block:: bash

   2017-08-26T10:35:10: +=================================+
   2017-08-26T10:35:10: | Maximum likelihood optimisation |
   2017-08-26T10:35:10: +=================================+
   2017-08-26T10:35:12:  >Iteration   0: -logL=-177232.327, Lambda=1.0e-03
   2017-08-26T10:35:12:  >Iteration   1: -logL=-198931.081, Lambda=1.0e-03, delta=21698.754, step=1.0e+00, max(|grad|)=-22713.401202 [CutoffEnergy:4]
   2017-08-26T10:35:12:  >Iteration   2: -logL=-200263.519, Lambda=1.0e-04, delta=1332.438, step=1.0e+00, max(|grad|)=-67648.133552 [CutoffEnergy:4]
   2017-08-26T10:35:12:  >Iteration   3: -logL=-212350.035, Lambda=1.0e-05, delta=12086.515, step=1.0e+00, max(|grad|)=23036.151920 [Index2:6]
   2017-08-26T10:35:12:  >Iteration   4: -logL=-213951.781, Lambda=1.0e-06, delta=1601.746, step=1.0e+00, max(|grad|)=-8297.648633 [CutoffEnergy:4]
   2017-08-26T10:35:13:  >Iteration   5: -logL=-215341.362, Lambda=1.0e-07, delta=1389.581, step=1.0e+00, max(|grad|)=2747.034517 [Index2:6]
   2017-08-26T10:35:13:  >Iteration   6: -logL=-215490.429, Lambda=1.0e-08, delta=149.067, step=1.0e+00, max(|grad|)=697.301361 [Index2:6]
   2017-08-26T10:35:13:   Iteration   7: -logL=-215490.429, Lambda=1.0e-09, delta=-14.065, step=1.0e+00, max(|grad|)=-1205.683279 [Index1:3] (stalled)
   2017-08-26T10:35:13:   Iteration   8: -logL=-215490.429, Lambda=1.0e-08, delta=-14.064, step=1.0e+00, max(|grad|)=-1205.659884 [Index1:3] (stalled)
   2017-08-26T10:35:14:   Iteration   9: -logL=-215490.429, Lambda=1.0e-07, delta=-14.055, step=1.0e+00, max(|grad|)=-1205.425967 [Index1:3] (stalled)
   2017-08-26T10:35:14:   Iteration  10: -logL=-215490.429, Lambda=1.0e-06, delta=-13.963, step=1.0e+00, max(|grad|)=-1203.090620 [Index1:3] (stalled)
   2017-08-26T10:35:14:   Iteration  11: -logL=-215490.429, Lambda=1.0e-05, delta=-13.066, step=1.0e+00, max(|grad|)=-1180.113526 [Index1:3] (stalled)
   2017-08-26T10:35:14:   Iteration  12: -logL=-215490.429, Lambda=1.0e-04, delta=-6.263, step=1.0e+00, max(|grad|)=-983.079153 [Index1:3] (stalled)
   2017-08-26T10:35:14:  >Iteration  13: -logL=-215496.367, Lambda=1.0e-03, delta=5.937, step=1.0e+00, max(|grad|)=-285.734284 [Index1:3]
   2017-08-26T10:35:15:   Iteration  14: -logL=-215496.367, Lambda=1.0e-04, delta=-8.812, step=1.0e+00, max(|grad|)=-731.894749 [Index1:3] (stalled)
   2017-08-26T10:35:15:  >Iteration  15: -logL=-215497.536, Lambda=1.0e-03, delta=1.169, step=1.0e+00, max(|grad|)=-153.668138 [Index1:3]
   2017-08-26T10:35:15:   Iteration  16: -logL=-215496.136, Lambda=1.0e-04, delta=-1.400, step=1.0e+00, max(|grad|)=-328.522942 [Index1:3] (stalled)
   2017-08-26T10:35:15:  >Iteration  17: -logL=-215498.376, Lambda=1.0e-03, delta=2.240, step=1.0e+00, max(|grad|)=3.257091 [Index2:6]
   2017-08-26T10:35:16:  >Iteration  18: -logL=-215498.383, Lambda=1.0e-04, delta=0.007, step=1.0e+00, max(|grad|)=-8.270351 [CutoffEnergy:4]
   2017-08-26T10:35:16:  >Iteration  19: -logL=-215498.384, Lambda=1.0e-05, delta=0.002, step=1.0e+00, max(|grad|)=0.637485 [Index2:6]
   2017-08-26T10:35:16:
   2017-08-26T10:35:16: +=========================================+
   2017-08-26T10:35:16: | Maximum likelihood optimisation results |
   2017-08-26T10:35:16: +=========================================+
   2017-08-26T10:35:16: === GOptimizerLM ===
   2017-08-26T10:35:16:  Optimized function value ..: -215498.384
   2017-08-26T10:35:16:  Absolute precision ........: 0.005
   2017-08-26T10:35:16:  Acceptable value decrease .: 2
   2017-08-26T10:35:16:  Optimization status .......: converged
   2017-08-26T10:35:16:  Number of parameters ......: 14
   2017-08-26T10:35:16:  Number of free parameters .: 5
   2017-08-26T10:35:16:  Number of iterations ......: 19
   2017-08-26T10:35:16:  Lambda ....................: 1e-06
   2017-08-26T10:35:16:  Maximum log likelihood ....: 215498.384
   2017-08-26T10:35:16:  Observed events  (Nobs) ...: 202330.000
   2017-08-26T10:35:16:  Predicted events (Npred) ..: 202318.082 (Nobs - Npred = 11.9178772509331)
   2017-08-26T10:35:16: === GModels ===
   2017-08-26T10:35:16:  Number of models ..........: 3
   2017-08-26T10:35:16:  Number of parameters ......: 14
   2017-08-26T10:35:16: === GModelSky ===
   2017-08-26T10:35:16:  Name ......................: Vela
   2017-08-26T10:35:16:  Instruments ...............: all
   2017-08-26T10:35:16:  Instrument scale factors ..: unity
   2017-08-26T10:35:16:  Observation identifiers ...: all
   2017-08-26T10:35:16:  Model type ................: PointSource
   2017-08-26T10:35:16:  Model components ..........: "PointSource" * "SuperExponentialCutoffPowerLaw" * "Constant"
   2017-08-26T10:35:16:  Number of parameters ......: 8
   2017-08-26T10:35:16:  Number of spatial par's ...: 2
   2017-08-26T10:35:16:   RA .......................: 128.84 [-360,360] deg (fixed,scale=1)
   2017-08-26T10:35:16:   DEC ......................: -45.18 [-90,90] deg (fixed,scale=1)
   2017-08-26T10:35:16:  Number of spectral par's ..: 5
   2017-08-26T10:35:16:   Prefactor ................: 4.3553344038324e-09 +/- 4.75309654008435e-10 [1e-16,1e-06] ph/cm2/s/MeV (free,scale=1e-09,gradient)
   2017-08-26T10:35:16:   Index1 ...................: -1.3435116067011 +/- 0.0311270092486095 [-0,-5]  (free,scale=-1,gradient)
   2017-08-26T10:35:16:   CutoffEnergy .............: 988.499237722443 +/- 181.86709257978 [1,1000000] MeV (free,scale=1000,gradient)
   2017-08-26T10:35:16:   PivotEnergy ..............: 1000 [10,1000000] MeV (fixed,scale=1000,gradient)
   2017-08-26T10:35:16:   Index2 ...................: 0.587129723482927 +/- 0.0279303486200403 [0.1,5]  (free,scale=1,gradient)
   2017-08-26T10:35:16:  Number of temporal par's ..: 1
   2017-08-26T10:35:16:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2017-08-26T10:35:16: === GModelSky ===
   2017-08-26T10:35:16:  Name ......................: Galactic_diffuse
   2017-08-26T10:35:16:  Instruments ...............: all
   2017-08-26T10:35:16:  Instrument scale factors ..: unity
   2017-08-26T10:35:16:  Observation identifiers ...: all
   2017-08-26T10:35:16:  Model type ................: DiffuseSource
   2017-08-26T10:35:16:  Model components ..........: "DiffuseMapCube" * "Constant" * "Constant"
   2017-08-26T10:35:16:  Number of parameters ......: 3
   2017-08-26T10:35:16:  Number of spatial par's ...: 1
   2017-08-26T10:35:16:   Normalization ............: 1 [0.1,10]  (fixed,scale=1,gradient)
   2017-08-26T10:35:16:  Number of spectral par's ..: 1
   2017-08-26T10:35:16:   Normalization ............: 1.10915296009432 +/- 0.00503367301009883 [0.1,1000]  (free,scale=1,gradient)
   2017-08-26T10:35:16:  Number of temporal par's ..: 1
   2017-08-26T10:35:16:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2017-08-26T10:35:16: === GModelSky ===
   2017-08-26T10:35:16:  Name ......................: Extragalactic_diffuse
   2017-08-26T10:35:16:  Instruments ...............: all
   2017-08-26T10:35:16:  Instrument scale factors ..: unity
   2017-08-26T10:35:16:  Observation identifiers ...: all
   2017-08-26T10:35:16:  Model type ................: DiffuseSource
   2017-08-26T10:35:16:  Model components ..........: "DiffuseIsotropic" * "FileFunction" * "Constant"
   2017-08-26T10:35:16:  Number of parameters ......: 3
   2017-08-26T10:35:16:  Number of spatial par's ...: 1
   2017-08-26T10:35:16:   Value ....................: 1 [0,10]  (fixed,scale=1,gradient)
   2017-08-26T10:35:16:  Number of spectral par's ..: 1
   2017-08-26T10:35:16:   Normalization ............: 1 [0,1000]  (fixed,scale=1,gradient)
   2017-08-26T10:35:16:  Number of temporal par's ..: 1
   2017-08-26T10:35:16:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)


The fit converged after 19 iterations with spectral parameters that are
reasonably close to those found in
`Abdo et al. 2010, ApJ, 713, 154 <http://iopscience.iop.org/article/10.1088/0004-637X/713/1/154/pdf>`_
using the Fermi-LAT Science Tools.

  .. note::

     The fit results can be compared to the corresponding output of gtlike
     when the tool is applied to the same test data. As you will see below,
     the fitted model parameters are basically identical.

     .. code-block:: bash

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

     .. code-block:: bash

        2017-08-26T13:05:12: === GModelSky ===
        2017-08-26T13:05:12:  Name ......................: Vela
        2017-08-26T13:05:12:  Instruments ...............: all
        2017-08-26T13:05:12:  Test Statistic ............: 326326.556902542

     which is also close to the value determined by gtlike.
