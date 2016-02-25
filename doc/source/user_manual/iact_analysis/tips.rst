.. _sec_tips:

Tips, Tricks and FAQs
=====================
This section shows some of the wrinkles that could fascilitate the users' life when running analyses.

Visualise observations
----------------------
A very common task is to inspect the observation container and visualise several quantities.

Create a skymap of counts
^^^^^^^^^^^^^^^^^^^^^^^^^
The tool :ref:`ctskymap` allows to create a single map of counts contained inside an observation container.
In contrast to :ref:`ctbin`, ctskymap takes into account all events that fall into the skymap energy range and spatial range.
Since :ref:`ctbin` is indended for binned/stacked analysis, it respects energy thresholds, RoI selections etc.
Accordingly, :ref:`ctskymap` should be used for a simple visualisation of the counts.


.. code-block:: bash

	$ ctskymap debug=yes
	Input event list or observation definition XML file [events.fits] obs.xml 
	Lower energy limit (TeV) [0.1] 
	Upper energy limit (TeV) [100.0] 
	First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 
	Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 
	Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR] 
	Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] 
	Image scale (in degrees/pixel) [0.02] 
	Size of the X axis in pixels [200] 
	Size of the Y axis in pixels [200] 
	Output skymap file [skymap.fits] 
	
Create region file of pointings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In particular, one might want to know how the pointings of the observation container distribute across the sky.
For this purpose, :ref:`csobsinfo` provides the hidden parameter ``ds9file`` which signals to produce a region file that
contains these information:

.. code-block:: bash

	$ csobsinfo ds9file=pointings.reg
	Parfile "csobsinfo.par" not found. Create default parfile.
	Event list, counts cube, or observation definition file [obs.xml] 
	
Together with the previously produced skymap, one can plot the pointings using e.g. ds9:

.. code-block:: bash
  
  $ ds9 skymap.fits -regions pointings.reg

Plot zenith angle distribution of observations
----------------------------------------------
The tool :ref:`csobsinfo` dumps more information into the log file that might be of interest to catch in a script.
Therefore, one could simply use this tool via python and plot some distributions:

.. code-block:: python

  # Run script
  import cscripts
  info = cscripts.csobsinfo()
  info["inobs"] = "selected_obs.xml"
  info.run()
  
  # Plot data
  import matplotlib.pyplot as plt
  plt.hist(info.zeniths(), bins=10, range=[0,90])
  plt.show()
  
This will plot a histogram showing the zenith angle distribution of the observation container. In addition there are more methods
to access information in python:

.. code-block:: python

  # csobsinfo can also return:
  info.azimuths() # list of azimuth values
  info.offsets()  # list of offset values (only computed if hidden parameter offset=yes was specified)
  info.ebounds()  # gammalib.GEbounds object of energy ranges
  info.gti()      # gammalib.GGti object containing good time intervals 
  

Visualise models
----------------
Similar to :ref:`csobsinfo`, there is :ref:`csmodelinfo` to find out more about the content of a model XML file. In order to show the
position and sizes of the model on top of a skymap, this tool has the hidden parameter ``ds9file``, too.

.. code-block:: bash

	csmodelinfo ds9file=models.reg
	Input model XML file [$CTOOLS/share/models/crab.xml]

There are several options regarding the color, text and other attributes of the region file. To see a full list, 
:ref:`visit the reference page <csmodelinfo>`.
Analogous to the pointings of the observation container, the models can be visualised using e.g. ds9:

.. code-block:: bash
  
  $ ds9 skymap.fits -regions models.reg

Manipulating models in python
-----------------------------
This step will give some guidance how to work with model XML files in python and how to manipulate their content

.. code-block:: python

	import gammalib
	
	# Open model file
	models = gammalib.GModels("$CTOOLS/share/models/crab.xml")
	
	# Access a model component
	src = models["Crab"]
	print(src)
	
	# Retrieve the spectral or spatial component
	print(src.spectral())
	print(src.spatial())
	
	# Print the spectral parameter "Prefactor"
	prefactor = src.spectral()["Prefactor"]
	print("Prefactor:",prefactor.value(),"+-",prefactor.error())
	
	# Set the prefactor value
	prefactor.value(3.5e-17)
	
	# Loop over models and fix all parameters
	for model in models:
		
		# Loop over model parameters
		for par in model:
			
			# Fix parameter
			par.fix()
	
	# Release one specific parameter
	models["Crab"]["Prefactor"].free()
	
	# Set the parameter range
	models["Crab"]["Prefactor"].min(1e-18)
	models["Crab"]["Prefactor"].max(1e-16)
	
	# ... or in one step:
	models["Crab"]["Prefactor].range(1e-18, 1e-16)
	
	# Remove model from container
	models.remove("Crab")
	
	# Save model to another XML file
	models.save("mymodels.xml")
  

Retrieve likelihood values from :ref:`ctlike`
---------------------------------------------
For some purposes, it might be useful to retrieve the fitted likelihood value and other results of the fit with :ref:`ctlike`.

.. code-block:: python

	import gammalib
	import ctools
	
	# Create ctlike
	like = ctools.ctlike()
	like["inobs"]   = "selected_obs.xml"
	like["inmodel"] = "mymodels.xml"
	like.run()
	
	# Get optimizer from ctlike
	opt = like.opt()
	print(opt)
	
	# Get likelihood value from optimizer
	print(opt.value())
	
	# Get fit status:
	print(opt.status())
	# 0: converged
	# 1: stalled
	# 2: singular curvature matrix encountered
	# 3: curvature matrix not positive definite
	# 4: errors are inaccurate

Speed up analysis
-----------------
In some cases it may occur that the analysis takes very long. There are several reasons that can slow down the analysis quite severly.
Not all issues can be tackled. Nevertheless, here is a list of actions that have proven to speed up the fit.

* Consider switching to binned analysis if the observation time exceeds 50-100 hours.
* Diffuse and extended sources in the RoI are computational expensive. However, in unbinned analysis they impact the speed much more.
* Usage of energy dispersion in unbinned analysis could also cause a long fit. Try switching it off if speed is important.
* Switch to a machine that supports OpenMP. This allows to parallelise the fit onto several cores.

Compute excess maps
-------------------
A very important means to visualise images of sources are excess maps. The tool :ref:`csresmap` can also be used to create such a map.
The user, however, has to manually remove the source(s) of interest from the input model XML file.

.. code-block:: python

  import gammalib
  models = gammalib.GModels("crab_models.xml")
  models.remove("Crab")
  models.save("models_without_crab.xml")
  
Of course, this can also be done by editing the XML file with an editor.

Subsequently, :ref:`csresmap` can be executed using the ``algorithm=SUB`` parameter:

.. code-block:: bash

	$ csresmap
	Input event list, counts cube, or observation definition XML file [selected_obs.xml]
	Input model XML file [models_without_crab.xml]
	First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 
	Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 
	Coordinate System (CEL|GAL) [CEL] 
	Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR] 
	Size of the X axis in pixels [100]
	Size of the Y axis in pixels [100]
	Pixel size (deg/pixel) [0.02]
	Residual map computation algorithm (SUB|SUBDIV|SUBDIVSQRT) [SUB]
	Output residual map file [excessmap.fits] 


Compute flux maps
-----------------
The creation of a flux map is much more expensive in terms of computing than an excess map. In the section about
:ref:`high level analysis tools <sec_high_level>`, the usage of :ref:`cttsmap` was demonstrated. This tool can also
be used to generate flux maps. For this, the user has to incoorporate a test source in the XML model file. This source
should consist of a spectral component using a PowerLaw2 model:

.. code-block:: xml

	<source name="TestSource" type="PointSource">
		<spectrum type="PowerLaw2">
			<parameter scale="1e-07" name="Integral"   min="1e-07" max="1000.0"    value="1.0" free="1"/>
			<parameter scale="1.0"   name="Index"      min="-5.0"  max="+5.0"      value="-2.0" free="1"/>
			<parameter scale="1.0"   name="LowerLimit" min="10.0"  max="1000000.0" value="100.0" free="0"/>
			<parameter scale="1.0"   name="UpperLimit" min="10.0"  max="1000000.0" value="500000.0" free="0"/>
		</spectrum>
		<spatialModel type="SkyDirFunction">
			<parameter free="0" max="360" min="-360" name="RA" scale="1" value="83.6331" />
			<parameter free="0" max="90" min="-90" name="DEC" scale="1" value="22.0145" />
		</spatialModel>
	</source>

The further content of the rest of the XML model file depends on the user requirements:

* For an absolute flux map, all other sky models should be removed
* For a residual flux map, all other sky models should be kept

The tool :ref:`cttsmap` will create one skymap per free spectral parameter in the model. In the result FITS file,
there will be an extension called "Integral" that contains the flux map. For the above XML example, the source name
"TestSource" should be specified to :ref:`cttsmap`.

Speed up TS map computation
---------------------------
Since the computation of a TS map can be extremely time consuming, the option of splitting the computation into several jobs
is supported. This might be of particular interest if the user has access to a batch farm.

Split TS map computation
^^^^^^^^^^^^^^^^^^^^^^^^
For the purpose of job splitting, the hidden parameter ``binmin`` and ``binmax`` were included in the tool. These are integer parameter
that specify which bins should be computed. For instance, if the map should consist of 30x30 =90 pixels, the use could e.g. execute

.. code-block:: bash

	$ cttsmap binmin=0 binmax=299
	...
	$ cttsmap binmin=300 binmax=599
	...
	$ cttsmap binmin=600 binmax=899
	...
	
on different machines.

.. note::

  The output file name of individual jobs can be the same. If the parameters are specified, the output file name will be appended by the suffix
  "_<binmin>_<binmax>" to distinguish between the sliced TS maps.
  
Merge sliced TS maps
^^^^^^^^^^^^^^^^^^^^
The tool :ref:`cstsmapmerge` is intended to take care of merging the files that were produced while splitting the TS map computation.
There are several options to pass the files to be merged as arguments:

* a space-separated list of file names (e.g. tsmap1.fits tsmap2.fits)
* a comma-separated list of file names (e.g. tsmap1.fits,tsmap2.fits)
* a wildcard string (e.g. tsmap*.xml)
* an ASCII file containing the file names, one file per line (e.g. @mymaps.txt)

In this example the ASCII file method is presented:

.. code-block:: bash 

	# Put slice files into an ascii file
	$ ls tsmap_*.fits > tsmapfiles.txt
	
	# Run cstsmapmerge
	$ cstsmapmerge
	Input TS map FITS files [@tsmapfiles.txt]
	Output TS map FITS file [mytsmap.fits]


Creating a python analysis pipeline
-----------------------------------