.. _howto_comptel_fitting:

Fit a model to the data
-----------------------

  .. admonition:: What you will learn

     You will learn how to use :ref:`ctlike` to fit a parametric model to
     COMPTEL data.

Start with creating an
:ref:`observation definition file <glossary_obsdef>`
that defines the names of the counts cube files, geometry factor files
and exposure map files. In addition you will need background model files,
but since those are not distributed with the data you will use the geometry
factors also as background files, which is a reasonable zero order
approximation of the instrumental background distribution in COMPTEL.

Below the
:ref:`observation definition file <glossary_obsdef>`
that is needed for the following analyses. Since the four counts cubes
correspond to the four energy bands, and every observation holds a single
counts cube, an observation has been added for every energy band.
The ``IAQ`` parameter indicates the instrument response, which needs
to be given for each band in the specified format.

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <observation_list title="observation library">
     <observation name="Crab" id="100001" instrument="COM">
       <parameter name="DRE" file="m50438_dre.fits"/>
       <parameter name="DRB" file="m34997_drg.fits"/>
       <parameter name="DRG" file="m34997_drg.fits"/>
       <parameter name="DRX" file="m32171_drx.fits"/>
       <parameter name="IAQ" value="ENERG(0.75-1.0)MeV"/>
     </observation>
     <observation name="Crab" id="200001" instrument="COM">
       <parameter name="DRE" file="m50439_dre.fits"/>
       <parameter name="DRB" file="m34997_drg.fits"/>
       <parameter name="DRG" file="m34997_drg.fits"/>
       <parameter name="DRX" file="m32171_drx.fits"/>
       <parameter name="IAQ" value="ENERG(1.0-3.0)MeV"/>
     </observation>
     <observation name="Crab" id="300001" instrument="COM">
       <parameter name="DRE" file="m50440_dre.fits"/>
       <parameter name="DRB" file="m34997_drg.fits"/>
       <parameter name="DRG" file="m34997_drg.fits"/>
       <parameter name="DRX" file="m32171_drx.fits"/>
       <parameter name="IAQ" value="ENERG(3.0-10.0)MeV"/>
     </observation>
     <observation name="Crab" id="400001" instrument="COM">
       <parameter name="DRE" file="m50441_dre.fits"/>
       <parameter name="DRB" file="m34997_drg.fits"/>
       <parameter name="DRG" file="m34997_drg.fits"/>
       <parameter name="DRX" file="m32171_drx.fits"/>
       <parameter name="IAQ" value="ENERG(10.0-30.0)MeV"/>
     </observation>
   </observation_list>

Then you also need a
:ref:`model definition file <glossary_moddef>`
that describes the source and background model that you want to fit. In the
example below, a point source with a ower law is fit at the position of the
Crab on top of an instrumental background model.

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <source_library title="source library">
     <source name="Crab" type="PointSource" tscalc="1">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1e-3" value="1.0"  min="-1000.0" max="1000.0" free="1"/>
         <parameter name="Index"       scale="-1"   value="2.0"  min="-5.0"    max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1.0"  value="1.65" min="0.01"    max="1000.0" free="0"/>
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="0"/>
         <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="0"/>
       </spatialModel>
     </source>
     <source name="Background(0.75-1.0)MeV" type="DRBFitting" instrument="COM" id="100001">
       <node>
         <parameter name="Phibar"        scale="1.0" value="1.0"  min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
       </node>
       <node>
         <parameter name="Phibar"        scale="1.0" value="3.0"  min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
       </node>
       <node>
         <parameter name="Phibar"        scale="1.0" value="5.0"  min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
       </node>
       <node>
         <parameter name="Phibar"        scale="1.0" value="7.0"  min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
       </node>
       <node>
         <parameter name="Phibar"        scale="1.0" value="9.0"  min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
       </node>
       <node>
         <parameter name="Phibar"        scale="1.0" value="11.0" min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
       </node>
       <node>
         <parameter name="Phibar"        scale="1.0" value="13.0" min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
       </node>
       <node>
         <parameter name="Phibar"        scale="1.0" value="15.0" min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
       </node>
       <node>
         <parameter name="Phibar"        scale="1.0" value="17.0" min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="1.0"  min="0.0" max="1000.0" free="1"/>
       </node>
       ...
     </source>
     <source name="Background(1.0-3.0)MeV" type="DRBFitting" instrument="COM" id="200001">
       <node>
         <parameter name="Phibar"        scale="1.0" value="1.0"  min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
       </node>
       ...
     </source>
     <source name="Background(3.0-10.0)MeV" type="DRBFitting" instrument="COM" id="300001">
       <node>
         <parameter name="Phibar"        scale="1.0" value="1.0"  min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="1.0"  min="0.0" max="1000.0" free="1"/>
       </node>
       ...
     </source>
     <source name="Background(10.0-30.0)MeV" type="DRBFitting" instrument="COM" id="400001">
       <node>
         <parameter name="Phibar"        scale="1.0" value="1.0"  min="0.0" max="50.0"   free="0"/>
         <parameter name="Normalization" scale="1.0" value="1.0"  min="0.0" max="1000.0" free="1"/>
       </node>
       ...
     </source>
   </source_library>

Since the file is relatively long it is not entirely reproduced here, but you
can download the file :download:`models.xml <models.xml>` by clicking on
the file name.

The first source in the
:ref:`model definition file <glossary_moddef>`
is the point source for the Crab, the next four sources are background
model components for each of the energy bands.

You may have recognised the ``id`` attribute for each model component. This
attribute links a given model component to a given observation. For example,
``Background(1.0-3.0)MeV`` has ``id="200001"`` and will hence only be used
to model the background for the observation that corresponds to the 1 - 3 MeV
band.

Each background model implements a so called Phibar-fitting of the file that
was specified for the background model (recall that you specified the geometry
factors as background model). Phibar-fitting means that the background model
is fitted to the data for each Phibar-layer of the counts cube. Each background
model is composed of 25 nodes, corresponding to the number of Phibar layers
of a COMPTEL counts cube, and each node is defined by it's Phibar value
(a fixed parameter) and it's Normalization (in general a free parameter).
For some background models some normalizations of Phibar layers have been set
to zero and the parameters fixed to exclude layers from the fitting that do
not contain any data (the filling-factor of a counts cube is energy dependent;
for the 0.75 - 1 MeV energy bin the first 8 layers are in fact empty).

Now you are ready to fit the model using the :ref:`ctlike` model tool. You do
this as follows:

.. code-block:: bash

   $ ctlike
   Input event list, counts cube or observation definition XML file [events.fits] obs.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models.xml
   Output model definition XML file [crab_results.xml]

:ref:`ctlike` creates on output the
:ref:`model definition file <glossary_moddef>`
``crab_results.xml`` where the parameter values were updated by their
fitted values, and where the statistical uncertainties were added using the
``error`` attribute. To investigate how the fit went you should inspect the
log file ``ctlike.log`` that was also created by :ref:`ctlike`:

.. code-block:: bash

   2017-08-26T22:05:05: +=================================+
   2017-08-26T22:05:05: | Maximum likelihood optimisation |
   2017-08-26T22:05:05: +=================================+
   2017-08-26T22:05:06:  >Iteration   0: -logL=988086.091, Lambda=1.0e-03
   2017-08-26T22:05:07:  >Iteration   1: -logL=670534.792, Lambda=1.0e-03, delta=317551.299, step=1.0e+00, max(|grad|)=86219.859480 [Index:3]
   2017-08-26T22:05:08:  >Iteration   2: -logL=387046.862, Lambda=1.0e-04, delta=283487.930, step=1.0e+00, max(|grad|)=69886.869623 [Index:3]
   2017-08-26T22:05:08:  >Iteration   3: -logL=154426.362, Lambda=1.0e-05, delta=232620.499, step=1.0e+00, max(|grad|)=47994.333526 [Index:3]
   2017-08-26T22:05:09:  >Iteration   4: -logL=-12393.885, Lambda=1.0e-06, delta=166820.247, step=1.0e+00, max(|grad|)=24763.811238 [Index:3]
   2017-08-26T22:05:10:  >Iteration   5: -logL=-108271.655, Lambda=1.0e-07, delta=95877.770, step=1.0e+00, max(|grad|)=7544.090872 [Index:3]
   2017-08-26T22:05:11:  >Iteration   6: -logL=-145906.957, Lambda=1.0e-08, delta=37635.302, step=1.0e+00, max(|grad|)=1167.483336 [Index:3]
   2017-08-26T22:05:12:  >Iteration   7: -logL=-152898.212, Lambda=1.0e-09, delta=6991.255, step=1.0e+00, max(|grad|)=110.879761 [Index:3]
   2017-08-26T22:05:12:  >Iteration   8: -logL=-153217.134, Lambda=1.0e-10, delta=318.921, step=1.0e+00, max(|grad|)=-2.969681 [Scale factor 2:111]
   2017-08-26T22:05:13:  >Iteration   9: -logL=-153218.344, Lambda=1.0e-11, delta=1.210, step=1.0e+00, max(|grad|)=-0.047390 [Index:3]
   2017-08-26T22:05:14:  >Iteration  10: -logL=-153218.344, Lambda=1.0e-12, delta=0.000, step=1.0e+00, max(|grad|)=-0.001891 [Index:3]
   2017-08-26T22:05:15:
   2017-08-26T22:05:15: +====================================+
   2017-08-26T22:05:15: | Maximum likelihood re-optimisation |
   2017-08-26T22:05:15: +====================================+
   2017-08-26T22:05:15:  >Iteration   0: -logL=-152378.830, Lambda=1.0e-03
   2017-08-26T22:05:16:  >Iteration   1: -logL=-152759.255, Lambda=1.0e-03, delta=380.425, step=1.0e+00, max(|grad|)=-2.480864 [Scale factor 10:21]
   2017-08-26T22:05:16:  >Iteration   2: -logL=-152760.518, Lambda=1.0e-04, delta=1.264, step=1.0e+00, max(|grad|)=-0.019031 [Scale factor 10:21]
   2017-08-26T22:05:16:  >Iteration   3: -logL=-152760.518, Lambda=1.0e-05, delta=0.000, step=1.0e+00, max(|grad|)=-0.000001 [Scale factor 10:21]
   2017-08-26T22:05:16:
   2017-08-26T22:05:16: +============================================+
   2017-08-26T22:05:16: | Maximum likelihood re-optimisation results |
   2017-08-26T22:05:16: +============================================+
   2017-08-26T22:05:16: === GOptimizerLM ===
   2017-08-26T22:05:16:  Optimized function value ..: -152760.518
   2017-08-26T22:05:16:  Absolute precision ........: 0.005
   2017-08-26T22:05:16:  Acceptable value decrease .: 2
   2017-08-26T22:05:16:  Optimization status .......: converged
   2017-08-26T22:05:16:  Number of parameters ......: 200
   2017-08-26T22:05:16:  Number of free parameters .: 87
   2017-08-26T22:05:16:  Number of iterations ......: 3
   2017-08-26T22:05:16:  Lambda ....................: 1e-06
   2017-08-26T22:05:16:
   2017-08-26T22:05:16: +=========================================+
   2017-08-26T22:05:16: | Maximum likelihood optimisation results |
   2017-08-26T22:05:16: +=========================================+
   2017-08-26T22:05:16: === GOptimizerLM ===
   2017-08-26T22:05:16:  Optimized function value ..: -153218.344
   2017-08-26T22:05:16:  Absolute precision ........: 0.005
   2017-08-26T22:05:16:  Acceptable value decrease .: 2
   2017-08-26T22:05:16:  Optimization status .......: converged
   2017-08-26T22:05:16:  Number of parameters ......: 206
   2017-08-26T22:05:16:  Number of free parameters .: 89
   2017-08-26T22:05:16:  Number of iterations ......: 10
   2017-08-26T22:05:16:  Lambda ....................: 1e-13
   2017-08-26T22:05:16:  Maximum log likelihood ....: 153218.344
   2017-08-26T22:05:16:  Observed events  (Nobs) ...: 527595.000
   2017-08-26T22:05:16:  Predicted events (Npred) ..: 527594.000 (Nobs - Npred = 1.00004221417475)
   2017-08-26T22:05:16: === GModels ===
   2017-08-26T22:05:16:  Number of models ..........: 5
   2017-08-26T22:05:16:  Number of parameters ......: 206
   2017-08-26T22:05:16: === GModelSky ===
   2017-08-26T22:05:16:  Name ......................: Crab
   2017-08-26T22:05:16:  Instruments ...............: all
   2017-08-26T22:05:16:  Test Statistic ............: 915.650908945536
   2017-08-26T22:05:16:  Instrument scale factors ..: unity
   2017-08-26T22:05:16:  Observation identifiers ...: all
   2017-08-26T22:05:16:  Model type ................: PointSource
   2017-08-26T22:05:16:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
   2017-08-26T22:05:16:  Number of parameters ......: 6
   2017-08-26T22:05:16:  Number of spatial par's ...: 2
   2017-08-26T22:05:16:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
   2017-08-26T22:05:16:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
   2017-08-26T22:05:16:  Number of spectral par's ..: 3
   2017-08-26T22:05:16:   Prefactor ................: 0.00237517119227384 +/- 0.000107673211964802 [-1,1] ph/cm2/s/MeV (free,scale=0.001,gradient)
   2017-08-26T22:05:16:   Index ....................: -2.47153532398462 +/- 0.0500838601327289 [5,-5]  (free,scale=-1,gradient)
   2017-08-26T22:05:16:   PivotEnergy ..............: 1 [0.01,1000] MeV (fixed,scale=1,gradient)
   2017-08-26T22:05:16:  Number of temporal par's ..: 1
   2017-08-26T22:05:16:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

The fit converged after 10 iterations and did a second fit without the
source component to determine the Test Statistic value of the Crab.
The source is detected with a Test Statistic of 915.7 which corresponds to
a detection significance of about 30 sigma.

.. note::
   `Kuiper et al. 2001, A&A, 378, 918 <http://cdsads.u-strasbg.fr/abs/2001A%26A...378..918K>`_
   presented the most comprehensive analysis of the COMPTEL Crab observations
   and summarize in their Table 3 the gamma-ray intensity as function of
   energy for the unpulsed and pulsed component of source. By adding the
   intensities and fitting the resulting values we obtained a prefactor
   of 2.5e-3, which is very close to the prefactor of 2.4e-3 obtained
   in the analysis above. The spectral index of the Kuiper's data is
   -2.2, a bit harder than the index of -2.5 obtained in the fit.

.. note::
   The position of a point source can also be adjusted by :ref:`ctlike` for
   an analysis of COMPTEL data, hence you may set the corresponding parameters
   to ``free="1"`` to determine also the best fitting source position.

.. warning::
   ctools supports for the moment only the fitting of point sources for
   COMPTEL data. Other spatial shapes will be implemented in the future.
