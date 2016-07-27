.. _sec_simulation:

Simulating IACT data
====================
This section describes how to run simulations with current IACT data. In particular, a guide how to simulate data for a given set of IRFs, or a given set of
observations is provided. In addition, tools to compute sensitivity curves or pull distributions are presented.
All tools are based on the event list simulator :ref:`ctobssim`.

Simulations for a specific observation container 
------------------------------------------------
On default, :ref:`ctobssim` builds a new observation container by querying relevant user parameters. In case
an observation XML container is already in place, the hidden parameter ``inobs`` can be specified (:ref:`see here how to create an
observation container <sec_analysis>`). It might be useful to run :ref:`ctselect` beforehand to ensure that the simulation knows about
energy thresholds, field of view, etc. If an unselected observation container is used, the simulation tool will query for global energy
boundaries and a global RoI radius. The following example illustrates the usage of an already selected observation container.

.. code-block:: bash
 
	$ ctobssim inobs=selected_obs.xml
	Input model XML file [crab_models.xml] 
	Output event data file or observation definition XML file [sim_obs.xml] 
	
Note that :ref:`ctobssim` writes a new event list file per observation into the current directory. If the observation container
is large, it might be convenient to create a new folder and specify the ``prefix`` keyword (similar to :ref:`ctselect`):

.. code-block:: bash
  
  $ mkdir simulated
  $ ctobssim inobs=selected_obs.xml prefix=simulated/sim_events_
  ...
  
The newly created observation container can be used in the same way as a normal observation container, i.e. it can be analysed with
:ref:`ctlike` etc. 

.. note::

	On default, :ref:`ctobssim` uses the same seed in every execution, i.e. the tool will always return the same simulated event list, regardless
	when, where and how it was executed. To steer this the user can provide the hidden integer parameter ``seed``:
	
	.. code-block:: bash
	
		$ ctobssim inobs=selected_obs.xml seed=2
  
  This might be important if simulations are run in parallel.

Simulating observation time for a given set of IRFs
---------------------------------------------------
There are other use cases to run simulations, for instance when there is no observation container available yet. This might be the case e.g
for an observation proposal: one would want to use a certain set of IRFs to simulate e.g. 50 hours of data to check whether a source
of interest can be detected. For this purpose it is necessary to first find a set of suitable IRFs to be used for
simulations. This can be easily done using :ref:`csfindobs` and using the ``expression`` parameter. The following example extracts
observations with at a certain zenith angle range (20-25 degrees) and azimuth angle range (90-100 degrees):

.. code-block:: bash
	
	$ csfindobs expression="AZ_PNT>90.0&&AZ_PNT<100.0&&ZEN_PNT>20.0&&ZEN_PNT<25.0"
	Name of FITS production (Run csiactdata to view your options) [fits-prod-name]
	Right ascension [83.6331] NONE
	Runlist outfile [zen20_az90.lis]
	
Of course, further details could be specified, e.g. one might want to check for a certain range of muon efficieny, or any other parameter which impacts IRFs.
The output runlist subsequently has to be converted to an observation XML file using :ref:`csiactobs` (:ref:`see here how to execute this tool <sec_analysis>`).
In the following it is assumed that the observation XML file of the selected representative observations is produced and called ``myobs.xml``.

Create caldb entry
^^^^^^^^^^^^^^^^^^
Now, one of the observations in this container have to be added to the local calibration database. This database is used
by ctools to extract IRFs for simulations.

.. code-block:: bash

	$ csobs2caldb
	Input observation definition file [myobs.xml] 
	Response output name (e.g. Zenith50) [zen20_az90] 
	
Note that the response name has been specified to "zen20_az90". Any name can be chosen which represents the observation which was inserted into the database.
On default, :ref:`csobs2caldb` uses the first observation in the provided XML file. This can be modifyied by specifying the hidden integer parameter index, e.g.:

.. code-block:: bash

  $ csobs2caldb index=1
  
This would use the second entry instead of the first from the observation container to extract the IRFs.

.. note::

  The ``index`` parameter starts counting from zero.  

Check caldb entries
^^^^^^^^^^^^^^^^^^^
In order to check if the calibration database has successfully been added, the tool :ref:`cscaldb` can be used to inspect the calibration database:

.. code-block:: bash

	$ cscaldb debug=yes
	...
	2016-02-25T15:57:11: +==============+
	2016-02-25T15:57:11: | Mission: cta |
	2016-02-25T15:57:11: +==============+
	2016-02-25T15:57:11: === Response functions in database "hess" ===
	2016-02-25T15:57:11: zen20_az90
	2016-02-25T15:57:11: 
	2016-02-25T15:57:11: === Response functions in database "prod2" ===
	2016-02-25T15:57:11: North_0.5h
	2016-02-25T15:57:11: North_50h
	2016-02-25T15:57:11: North_5h
	2016-02-25T15:57:11: South_0.5h
	2016-02-25T15:57:11: South_50h
	2016-02-25T15:57:11: South_5h
	...

Note that the observation of the actual instrument ("hess" in this example) was added to the CTA mission. This is mainly for structural and simplicity reasons
and might be improved eventually. Nevertheless, this caldb entry can now be used to simulate any observation duration.

Simulate
^^^^^^^^

.. code-block:: bash
  
	$ ctobssim
	RA of pointing (degrees) (0-360) [83.63] 
	Dec of pointing (degrees) (-90-90) [22.01] 
	Radius of FOV (degrees) (0-180) [2.5] 
	Start time (MET in s) [0.0] 
	End time (MET in s) [18000.0] 
	Lower energy limit (TeV) [0.5] 
	Upper energy limit (TeV) [50] 
	Calibration database [hess] 
	Instrument response function [zen20_az90] 
	Input model XML file [input_model.xml] 
	Output event data file or observation definition XML file [simulated_events.fits]
	
Here, a total duration of 5 hours at the crab position was simulated using the newly created caldb entry. Note that the choice of the energy range is up to the user.
One should bear in mind that simulations only make sense in energy ranges where the IRFs are well-defined. The simulation
assumes one observation of defined duration here. Accordingly, only one FITS event list is written as output. The input XML model
must now contain only one background model and (of course) all sky models that are needed. Here is an example of an input model
e.g. for the Crab Nebula:

.. code-block:: xml

	<?xml version="1.0" standalone="no"?>
	<source_library title="source library">
	  <source name="Crab" type="PointSource">
	    <spectrum type="PowerLaw">
	       <parameter name="Prefactor"   scale="1e-17" value="3.5"  min="1e-07" max="1000.0" free="1"/>
	       <parameter name="Index"       scale="-1"    value="2.5"  min="0.0"   max="+5.0"   free="1"/>
	       <parameter name="PivotEnergy" scale="1e6"   value="1.0"  min="0.01"  max="1000.0" free="0"/>
	    </spectrum>
	    <spatialModel type="PointSource">
	      <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="0"/>
	      <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="0"/>
	    </spatialModel>
	  </source>
	  <source name="CTABackgroundModel" type="CTAIrfBackground" instrument="CTA">
	    <spectrum type="PowerLaw">	
	      <parameter name="Prefactor"   scale="1.0"  value="1.0"  min="1e-3" max="1e+3"   free="1"/>	
	      <parameter name="Index"       scale="1.0"  value="0.0"  min="-5.0" max="+5.0"   free="1"/>
	      <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
	    </spectrum>
	  </source>	
	</source_library>

Note that since the caldb entry (:ref:`csobs2caldb`) has been produced under the mission "CTA", the event list is now considered a CTA observation (regardless of the input instrument).
Therefore, it is important to have "CTA" specified in the instrument attribute of the background model XML component. Of course the addition of any sky model can be
done without a problem.

The output of this kind of simulation is a single FITS event list file that can be used for further analysis.
For instance one could run :ref:`ctlike` in the following way:

.. code-block:: bash

	$ ctlike
	Input event list, counts cube or observation definition XML file [simulated_events.fits] 
	Calibration database [hess] 
	Instrument response function [zen20_az90] 
	Input model XML file [input_model.xml] 
	Output model XML file [simulated_fit_results.xml] 
	
Note that :ref:`ctlike` automatically queries for the IRF database in case no full observation container is provided. Any other tool taking
unbinned observations as input will behave the same way. Now all means to run simulations with current IACTs are at hand.
The following sections will describe tools for running higher level simulation analysis.

Create pull distributions
-------------------------
Pull distributions inspect the capability of the science tools to find source parameters in the data that were injected into the simulations.
The tool :ref:`cspull` will run several times the sequence of simulating and fitting the data. After each sequence, the pull value of the
fitted parameters is written into an ASCII file. The pull value ``p`` of parameter ``X`` is defined as follows:

.. math::
  p = \frac{X_{\rm fitted} - X_{\rm true}}{\Delta X_{\rm fitted}}
  
Where :math:`\Delta X_{\rm fitted}` is the derived error of the fitted parameter :math:`X_{\rm fitted}`. After several trials, the pull value should follow
a normal distribution. Knowing this, the tool may help to distinguish if a source is detectable and can be normally reconstructed.
The tool can either work on an input observation container or on the caldb entry that was created above. The number of trials is queried by the tool, too.

Example: observation container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

	$ cspull inobs=selected_obs.xml 
	Number of energy bins (0=unbinned) [0] 
	Input model XML file [crab_results.xml] 
	Output pull distribution file [pull.dat] 
	Number of trials [100] 

Note that the tool asks for the energy binning. It will use an unbinned analysis if zero bins are provided. Otherwise
it will query for binning parameters and will convert the simulated observations into a binned/stacked observation by running
:ref:`ctbin`, :ref:`ctexpcube`, :ref:`ctpsfcube` and :ref:`ctbkgcube` in order before fitting the simulated data with :ref:`ctlike`.

Example: caldb entry
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

	$ cspull
	RA of pointing (deg) (0-360) [83.6331] 
	Dec of pointing (deg) (-90-90) [22.0145] 
	Duration (in s) [18000.0] 
	Lower energy limit (TeV) [0.5]
	Upper energy limit (TeV) [50]
	Calibration database [hess] 
	Instrument response function [zen20_az90] 
	Number of energy bins (0=unbinned) [0]
	Input model XML file [input_models.xml] 
	Output pull distribution file [pull.dat] 
	Number of trials [100] 

Again, the caldb entry ``zen20_az90`` which was created above is used to simulate 5 hours of data. The data is simulated and analysed 100
times in a row and the output is written into the file ``pull.dat``.

Visualise pull distributions
----------------------------
The output ASCII file of :ref:`cspull` can be inspected by the user using own scripts. As a starting point, the ``example`` folder
contains a script that can be used to visualise the pull distribution histogram:

.. code-block:: bash

	$ python $CTOOLS/examples/show_pull_histogram.py pull.dat Pull_Crab_Prefactor 50

Run the script without arguments to see the usage. In case a wrong parameter name was provided, the tool will print available parameters on the screen.
The last integer parameter defines the binning of the histogram. In addition, the tool shows the standard normal distribution which should
match the pull distribution for a high number of trials if everything went well.

Create sensitivity curve
------------------------
A very important means to illustrate if a source can be detected with a certain significance is a sensitivity curve.
The tool :ref:`cssens` can be used to compute either an integral, or a differential sensitivity curve for a given observation container.
Analogous to :ref:`cspull`, the tool works either on an unbinned observation container or using a caldb entry. 
The script :ref:`cssens` computes the differential or integrated CTA sensitivity using maximum likelihood fitting of a test source.
The differential sensitivity is determined for a number of energy bins, the integral sensitivity is determined for a number
of energy thresholds. The test source is fitted to simulated data using ctlike to determine it’s detection significance as a
function of source flux. The source flux is then varied until the source significance achieves a given level, specified by the 
(hidden) significance parameter ``sigma``. To damp variations between individual Monte Carlo simulations, a sliding average is
applied in the significance computation (controlled by hidden ``num_avg parameter``). Note that this procedure might take a while to compute.
To switch from diffferential to integral sensitivity, the hidden parameter ``type="Integral"`` has to be specified. It should be noted that
the results are written to a file ``sensitivity.dat``. To change the output file name, the hidden parameter ``outfile=mysensitivity.dat``
can be provided.


Example: observation container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

	$ cssens inobs=selected_obs.xml
	Input model XML file [crab_results.xml] 
	Source name [Crab] 
	Lower energy limit (TeV) [0.5] 
	Upper energy limit (TeV) [50] 
	Number of energy bins for differential sensitivity computation [21] 

Example: caldb entry
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

	$ cssens 
	Calibration database [hess] 
	Instrument response function [zen20_az90] 
	Effective exposure time (s) [18000.0] 
	Radius of ROI (deg) [2.5] 
	Input model XML file [input_models.xml] 
	Source name [Crab] 
	Lower energy limit (TeV) [0.5] 
	Upper energy limit (TeV) [50] 
	Number of energy bins for differential sensitivity computation [21]
	
Note that the input model in the latter case only requires to contain one background model (see section above). Using this setup,
the user can specify the observation time for the sensitivity curve. In this example, again, 5 hours are used.

Visualise sensitivity curves
----------------------------
In order to show a sensitivity curve, there is also a script in the example folder:

.. code-block:: bash

	$ python $CTOOLS/examples/show_sensitivity.py sensitivity.dat
	
This will display a sensitivity curve for the applied settings.
