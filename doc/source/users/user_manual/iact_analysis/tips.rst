.. _sec_tips:

Tips, Tricks and FAQs
=====================
This section shows some of the wrinkles that could ease your life when
running IACT analyses.


Visualise observations
----------------------
A very common task is to inspect the observation container and visualise several
quantities.

Create a skymap of counts
^^^^^^^^^^^^^^^^^^^^^^^^^
The tool :ref:`ctskymap` allows to create a single map of counts contained
inside an observation container. In contrast to :ref:`ctbin`, ctskymap takes
into account all events that fall into the skymap energy range and spatial
range. Since :ref:`ctbin` is intended for stacked analysis, it respects energy
thresholds, RoI selections, etc. In contrast, :ref:`ctskymap` may be used for
a simple visualisation of the counts.

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
	
Alternatively, you may start from a counts cube and collapse the maps into one
single map via Python:

.. code-block:: python
 
   >>> import gammalib
   >>> map = gammalib.GSkyMap("myCountCube.fits") # Load counts cube
   >>> map.stack_maps()                           # Stack all cube maps into a single map
   >>> map.save("myCountMap.fits")                # Save map to file
	
Create region file of pointings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You may want to know how the pointings of the observation container distribute
across the sky. For this purpose, :ref:`csobsinfo` provides the hidden parameter
``ds9file`` which instructs the script produce a region file that contains the
pointing information:

.. code-block:: bash

   $ csobsinfo ds9file=pointings.reg debug=yes
   Event list, counts cube, or observation definition file [obs.xml]
	
Together with the previously produced skymap, you can plot the pointings using
for example ds9:

.. code-block:: bash
  
   $ ds9 skymap.fits -regions pointings.reg


Plot zenith angle distribution of observations
----------------------------------------------
The script :ref:`csobsinfo` logs information about an observation container
into a log file. You may access the information through dedicated class methods
and plot the information using for example matplotlib:

.. code-block:: python

   >>> # Run script
   >>> import cscripts
   >>> info = cscripts.csobsinfo()
   >>> info["inobs"] = "selected_obs.xml"
   >>> info.run()

   >>> # Plot data
   >>> import matplotlib.pyplot as plt
   >>> plt.hist(info.zeniths(), bins=10, range=[0,90])
   >>> plt.show()
  
This example will plot a histogram showing the zenith angle distribution of the
observation container. In addition there are more methods to access information
in Python:

.. code-block:: python

   >>> info.azimuths() # list of azimuth values
   >>> info.offsets()  # list of offset values (only computed if hidden parameter offset=yes was specified)
   >>> info.ebounds()  # gammalib.GEbounds object of energy ranges
   >>> info.gti()      # gammalib.GGti object containing good time intervals
  
There is also an example script to plot information about an observation
container which might be useful for a quick visualisation:

.. code-block:: bash
  
   $ python $CTOOLS/examples/show_obs.py selected_obs.xml


Visualise models
----------------
Similar to :ref:`csobsinfo`, there is the :ref:`csmodelinfo` script to find out
more about the content of a model XML file. In order to show the position and
sizes of the model on top of a skymap, this tool has the hidden parameter
``ds9file``, too.

.. code-block:: bash

   $ csmodelinfo ds9file=models.reg debug=yes
   Input model XML file [$CTOOLS/share/models/crab.xml]

There are several options regarding the color, text and other attributes of the
region file. To see a full list, :ref:`visit the reference page <csmodelinfo>`.
Analogous to the pointings of the observation container, the models can be
visualised using for example ds9:

.. code-block:: bash
  
   $ ds9 skymap.fits -regions models.reg


Manipulating models in python
-----------------------------
This example will give you some guidance on how to work with model XML files in
Python and how to manipulate their content:

.. code-block:: python

   >>> import gammalib

   >>> # Open model file
   >>> models = gammalib.GModels("$CTOOLS/share/models/crab.xml")

   >>> # Access a model component
   >>> src = models["Crab"]
   >>> print(src)

   >>> # Retrieve the spectral or spatial component
   >>> print(src.spectral())
   >>> print(src.spatial())

   >>> # Print the spectral parameter "Prefactor"
   >>> prefactor = src.spectral()["Prefactor"]
   >>> print("Prefactor: "+str(prefactor.value())+" +- "+str(prefactor.error()))

   >>> # Set the prefactor value
   >>> prefactor.value(3.5e-17)

   >>> # Loop over models and fix all parameters
   >>> for model in models:
   >>>     for par in model:   # Loop over model parameters
   >>>         par.fix()       # Fix parameter

   >>> # Release one specific parameter
   >>> models["Crab"]["Prefactor"].free()

   >>> # Set the parameter range
   >>> models["Crab"]["Prefactor"].min(1e-18)
   >>> models["Crab"]["Prefactor"].max(1e-16)

   >>> # ... or in one step:
   >>> models["Crab"]["Prefactor"].range(1e-18, 1e-16)

   >>> # Remove model from container
   >>> models.remove("Crab")

   >>> # Save model to another XML file
   >>> models.save("mymodels.xml")
  

Retrieve likelihood values from ctlike
--------------------------------------
For some purposes, it might be useful to retrieve the fitted likelihood value
and other results of the fit with :ref:`ctlike`.

.. code-block:: python

   >>> import gammalib
   >>> import ctools

   >>> # Create and run ctlike
   >>> like = ctools.ctlike()
   >>> like["inobs"]   = "selected_obs.xml"
   >>> like["inmodel"] = "mymodels.xml"
   >>> like.run()

   >>> # Get optimizer from ctlike
   >>> opt = like.opt()
   >>> print(opt)

   >>> # Get likelihood value from optimizer
   >>> print(opt.value())

   >>> # Get fit status:
   >>> print(opt.status())
   >>> # 0: converged
   >>> # 1: stalled
   >>> # 2: singular curvature matrix encountered
   >>> # 3: curvature matrix not positive definite
   >>> # 4: errors are inaccurate


Speed up analysis
-----------------
In some cases it may occur that the analysis takes very long. There are several
reasons that can slow down the analysis quite severly. Not all issues can be
tackled. Nevertheless, here is a list of actions that have proven to speed up
the fit.

* Consider switching to stacked analysis if the observation time exceeds 50-100
  hours.
* Diffuse and extended sources in the RoI are computational expensive. However,
  in unbinned analysis they impact the speed much more.
* Usage of energy dispersion could also cause a long fit. Try switching it off
  if speed is important. Keep in mind this will have an effect on the results.
* Switch to a machine that supports OpenMP. This allows to parallelise the fit
  onto several cores.


Compute excess maps
-------------------
A very important means to visualise images of sources are excess maps. The script
:ref:`csresmap` can also be used to create such a map. The user, however, has to
manually remove the source(s) of interest from the input model XML file.

.. code-block:: python

   >>> import gammalib
   >>> models = gammalib.GModels("crab_models.xml")
   >>> models.remove("Crab")
   >>> models.save("models_without_crab.xml")
  
Of course, this can also be done by editing the XML file with an editor.

Subsequently, :ref:`csresmap` can be executed using the ``algorithm=SUB``
parameter:

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
The creation of a flux map is much more expensive in terms of computing than an
excess map. In the section about :ref:`high level analysis tools <sec_high_level>`,
the usage of :ref:`cttsmap` was demonstrated. This tool can also be used to
generate flux maps. For this, the user has to incorporate a test source in the
XML model file. This source should consist of a spectral component using a
``PowerLawPhotonFlux`` model:

.. code-block:: xml

	<source name="TestSource" type="PointSource">
	 <spectrum type="PowerLawPhotonFlux">
	  <parameter scale="1e-07" name="PhotonFlux" min="1e-07" max="1000.0"    value="1.0" free="1"/>
	  <parameter scale="1.0"   name="Index"      min="-5.0"  max="+5.0"      value="-2.0" free="1"/>
	  <parameter scale="1.0"   name="LowerLimit" min="10.0"  max="1000000.0" value="100.0" free="0"/>
	  <parameter scale="1.0"   name="UpperLimit" min="10.0"  max="1000000.0" value="500000.0" free="0"/>
	 </spectrum>
	 <spatialModel type="PointSource">
	  <parameter free="0" max="360" min="-360" name="RA"  scale="1" value="83.6331" />
	  <parameter free="0" max="90"  min="-90"  name="DEC" scale="1" value="22.0145" />
	 </spatialModel>
	</source>

The further content of the rest of the XML model file depends on the user
requirements:

* For an absolute flux map, all other sky models should be removed
* For a residual flux map, all other sky models should be kept

The tool :ref:`cttsmap` will create one skymap per free spectral parameter in
the model. In the result FITS file, there will be an extension called ``Integral``
that contains the flux map. For the above XML example, the source name
``TestSource`` should be specified to :ref:`cttsmap`.


Speed up TS map computation
---------------------------
Since the computation of a TS map can be extremely time consuming, the option of
splitting the computation into several jobs is supported. This might be of
particular interest if the user has access to a batch farm.

Split TS map computation
^^^^^^^^^^^^^^^^^^^^^^^^
For the purpose of job splitting, the hidden parameter ``binmin`` and ``binmax``
were included in the tool. These are integer parameter that specify which bins
should be computed. For instance, if the map should consist of 30x30(=900) pixels,
the user could for example execute

.. code-block:: bash

	$ cttsmap binmin=0 binmax=299 outmap=tsmap_0_299.fits
	...
	$ cttsmap binmin=300 binmax=599 outmap=tsmap_300_599.fits
	...
	$ cttsmap binmin=600 binmax=899 outmap=tsmap_600_899.fits
	...
	
Each command could run on a different machine.

.. note::

  The output file name of individual jobs should be different. Otherwise files
  could overwrite each other. The naming of the individual slices is up to the
  user.
  
The script :ref:`cstsmapsplit` will take care of this bookkeeping. It creates an
ASCII file containing all the commands according to the user input. The following
example will create split the computation of a TS map into 2000 separate task.
Each task will only compute 5 bins on its own. This is very useful if the
observation container is large or the fit simply takes a long time.

.. code-block:: bash

  $ cstsmapsplit
  Input event list, counts cube or observation definition XML file [selected_obs.xml]
  Input model definition XML file [$CTOOLS/share/models/crab.xml]
  First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 
  Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 
  Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR] 
  Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] 
  Image scale (in degrees/pixel) [0.02]
  Size of the X axis in pixels [100]
  Size of the Y axis in pixels [100]
  Test source name [Crab] 
  Output Test Statistic map file [tsmap.fits] 
  Number of TS map bins per task [5] 
  Compute null hypothesis first? [yes] 
  ASCII file containing all commands [commands.dat] 
  
.. note::
  
  One can decide if we want to compute the null hypothesis first. This way, a
  parameter optimisation will be performed without the test source. The obtained
  likelihood value will then be passed to each individual task of :ref:`cttsmap`.
  If ``compute_null=no``, each task has to compute the null hypothesis itself.
  
Merge splitted TS maps
^^^^^^^^^^^^^^^^^^^^^^
The script :ref:`cstsmapmerge` is intended to take care of merging the files
that were produced while splitting the TS map computation. There are several
options to pass the files to be merged as arguments:

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
It is easily possible to build an own analysis workflow with a simple python
script. The following source code example shows a python script running from
gathering observations until fitting spectral points without storing intermediate
data products on disk. It assumes that the environment variable ``$VHEFITS`` is
set to the path where IACT FITS data is located.

.. code-block:: python
  
   >>> import gammalib
   >>> import ctools
   >>> import cscripts

   >>> # Set debug flag
   >>> debug = True

   >>> # Set flag to use energy dispersion
   >>> edisp = False

   >>> # Set inmodel file name
   >>> inmodel = "$GAMMALIB/test/data/model_point_plaw.xml"

   >>> # Expand environment variable
   >>> inmodel = gammalib.expand_env(inmodel)

   >>> # Extract coordinates and model properties
   >>> models  = gammalib.GModels(inmodel)
   >>> srcname = models[0].name()
   >>> ra      = models[0]["RA"].value()
   >>> dec     = models[0]["DEC"].value()

   >>> # Find FITS production name
   >>> iactdata          = cscripts.csiactdata()
   >>> iactdata["debug"] = debug
   >>> iactdata.run()

   >>> # Use first available production
   >>> prodname = iactdata.names()[0]

   >>> # Run csfindobs
   >>> findobs             = cscripts.csfindobs()
   >>> findobs["ra"]       = ra
   >>> findobs["dec"]      = dec
   >>> findobs["rad"]      = 2.5
   >>> findobs["prodname"] = prodname
   >>> findobs["debug"]    = debug
   >>> findobs["outfile"]  = "NONE"
   >>> findobs.run()

   >>> # Retrieve obervation IDs (runlist)
   >>> obs_ids = findobs.obs_ids()

   >>> # Build observation container
   >>> iactobs             = cscripts.csiactobs()
   >>> iactobs["prodname"] = prodname
   >>> iactobs["inmodel"]  = inmodel
   >>> iactobs["bkgpars"]  = 1
   >>> iactobs["outobs"]   = "NONE"
   >>> iactobs["outmodel"] = "NONE"
   >>> iactobs["debug"]    = debug
   >>> iactobs.runlist(obs_ids)
   >>> iactobs.run()

   >>> # Retrieve observation container and energy boundaries
   >>> obs     = iactobs.obs()
   >>> ebounds = iactobs.ebounds()

   >>> # Run ctselect
   >>> select             = ctools.ctselect(obs)
   >>> select["usepnt"]   = True
   >>> select["rad"]      = 2.5
   >>> select["usethres"] = "DEFAULT"
   >>> select["tmin"]     = 0.0
   >>> select["tmax"]     = 0.0
   >>> select["emin"]     = 0.1
   >>> select["emax"]     = 100.0
   >>> select["debug"]    = debug
   >>> select.run()

   >>> # Pass selected observations to ctlike
   >>> like          = ctools.ctlike(select.obs())
   >>> like["debug"] = debug
   >>> like["edisp"] = edisp
   >>> like.run()

   >>> # Compute a spectrum and save
   >>> spec             = cscripts.csspec(like.obs())
   >>> spec["srcname"]  = srcname
   >>> spec["emin"]     = ebounds.emin().TeV()
   >>> spec["emax"]     = ebounds.emax().TeV()
   >>> spec["enumbins"] = 10
   >>> spec["edisp"]    = edisp
   >>> spec["ebinalg"]  = "LOG"
   >>> spec["debug"]    = debug
   >>> spec["outfile"]  = "spectrum.fits"
   >>> spec.execute()

   >>> # Remove model for excess map computation
   >>> like.obs().models().remove(srcname)

   >>> # Compute an excess map
   >>> resmap             = cscripts.csresmap(like.obs())
   >>> resmap["xref"]     = ra
   >>> resmap["yref"]     = dec
   >>> resmap["proj"]     = "CAR"
   >>> resmap["coordsys"] = "CEL"
   >>> resmap["emin"]     = ebounds.emin().TeV()
   >>> resmap["emax"]     = ebounds.emax().TeV()
   >>> resmap["nxpix"]    = 100
   >>> resmap["nypix"]    = 100
   >>> resmap["binsz"]    = 0.02
   >>> resmap["outmap"]   = "excessmap.fits"
   >>> resmap.execute()

The results of this analysis workflow can be inspected using an example Python
script and for example ds9:

.. code-block:: bash
  
   $ python $CTOOLS/examples/show_spectrum.py spectrum.fits
   $ ds9 excessmap.fits
  