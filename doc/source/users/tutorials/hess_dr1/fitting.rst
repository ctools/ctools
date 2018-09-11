.. _hess_dr1_fitting:

Fitting the model components to the data
----------------------------------------

  .. admonition:: What you will learn

     You will learn **how to fit a parametric model to the data**.

Adjust source model name and set pivot energy to 1 TeV to bring it energy
interval of data.
Remove ``error`` attributes.

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Crab" type="PointSource">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   value="1" scale="5.7e-18" min="0" free="1" />
         <parameter name="Index"       value="1" scale="-2.48" min="-4" max="4" free="1" />
         <parameter name="PivotEnergy" value="1" scale="1000000" free="0" />
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA"  value="83.6192131308071" scale="1" free="1" />
         <parameter name="DEC" value="22.0199996472185" scale="1" free="1" />
       </spatialModel>
     </source>
     <source name="Background" type="CTAAeffBackground">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   value="1" scale="1e-13" min="0" free="1" />
         <parameter name="Index"       value="1" scale="-2.5"  min="-4" max="4" free="1" />
         <parameter name="PivotEnergy" value="1" scale="1000000" free="0" />
       </spectrum>
     </source>
   </source_library>


.. code-block:: bash

   $ ctlike debug=yes
   Input event list, counts cube or observation definition XML file [events.fits] obs_crab_selected.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml] crab_models.xml
   Output model definition XML file [crab_results.xml]

.. code-block:: none

   2018-09-11T21:17:31: +=========================================+
   2018-09-11T21:17:31: | Maximum likelihood optimisation results |
   2018-09-11T21:17:31: +=========================================+
   2018-09-11T21:17:31: === GOptimizerLM ===
   2018-09-11T21:17:31:  Optimized function value ..: 98422.688
   2018-09-11T21:17:31:  Absolute precision ........: 0.005
   2018-09-11T21:17:31:  Acceptable value decrease .: 2
   2018-09-11T21:17:31:  Optimization status .......: converged
   2018-09-11T21:17:31:  Number of parameters ......: 10
   2018-09-11T21:17:31:  Number of free parameters .: 6
   2018-09-11T21:17:31:  Number of iterations ......: 24
   2018-09-11T21:17:31:  Lambda ....................: 0.1
   2018-09-11T21:17:31:  Maximum log likelihood ....: -98422.688
   2018-09-11T21:17:31:  Observed events  (Nobs) ...: 9675.000
   2018-09-11T21:17:31:  Predicted events (Npred) ..: 9674.870 (Nobs - Npred = 0.129601074213497)
   2018-09-11T21:17:31: === GModels ===
   2018-09-11T21:17:31:  Number of models ..........: 2
   2018-09-11T21:17:31:  Number of parameters ......: 10
   2018-09-11T21:17:31: === GModelSky ===
   2018-09-11T21:17:31:  Name ......................: Crab
   2018-09-11T21:17:31:  Instruments ...............: all
   2018-09-11T21:17:31:  Instrument scale factors ..: unity
   2018-09-11T21:17:31:  Observation identifiers ...: all
   2018-09-11T21:17:31:  Model type ................: PointSource
   2018-09-11T21:17:31:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
   2018-09-11T21:17:31:  Number of parameters ......: 6
   2018-09-11T21:17:31:  Number of spatial par's ...: 2
   2018-09-11T21:17:31:   RA .......................: 83.6230374464172 +/- 0.00241074887155303 deg (free,scale=1)
   2018-09-11T21:17:31:   DEC ......................: 22.0249049416179 +/- 0.0022137104783551 deg (free,scale=1)
   2018-09-11T21:17:31:  Number of spectral par's ..: 3
   2018-09-11T21:17:31:   Prefactor ................: 4.81375364083447e-17 +/- 2.62970930530696e-18 [0,infty[ ph/cm2/s/MeV (free,scale=5.7e-18,gradient)
   2018-09-11T21:17:31:   Index ....................: -2.71079050853091 +/- 0.0652734279744327 [9.92,-9.92]  (free,scale=-2.48,gradient)
   2018-09-11T21:17:31:   PivotEnergy ..............: 1000000 MeV (fixed,scale=1000000,gradient)
   2018-09-11T21:17:31:  Number of temporal par's ..: 1
   2018-09-11T21:17:31:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2018-09-11T21:17:31: === GCTAModelAeffBackground ===
   2018-09-11T21:17:31:  Name ......................: Background
   2018-09-11T21:17:31:  Instruments ...............: all
   2018-09-11T21:17:31:  Instrument scale factors ..: unity
   2018-09-11T21:17:31:  Observation identifiers ...: all
   2018-09-11T21:17:31:  Model type ................: "PowerLaw" * "Constant"
   2018-09-11T21:17:31:  Number of parameters ......: 4
   2018-09-11T21:17:31:  Number of spectral par's ..: 3
   2018-09-11T21:17:31:   Prefactor ................: 1.80502240642969e-13 +/- 2.34991048262393e-15 [0,infty[ ph/cm2/s/MeV (free,scale=1e-13,gradient)
   2018-09-11T21:17:31:   Index ....................: -2.58303579414564 +/- 0.0128891124984818 [10,-10]  (free,scale=-2.5,gradient)
   2018-09-11T21:17:31:   PivotEnergy ..............: 1000000 MeV (fixed,scale=1000000,gradient)
   2018-09-11T21:17:31:  Number of temporal par's ..: 1
   2018-09-11T21:17:31:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
