.. _sec_high_level:

High level analysis tools for IACT data
=======================================

* :ref:`csspec` - calculate spectral points
* :ref:`cslightcrv` - calculate a light curve
* :ref:`csresmap` - compute a residual map
* :ref:`cttsmap` - compute a TS map (and flux map)

Inspecting observation definition files
---------------------------------------
After having created an observation definition file, the user might want to know a summary of its content.
Important quantities usually are the number of observations, the total livetime, or the mean zenith angle.
For this purpose, :ref:`csobsinfo` inspects an observation XML container for its content.

.. code-block:: bash

  $ csobsinfo
  Event list, counts cube, or observation definition file [obs.xml] selected_obs.xml 
  
This dumps general information of the observation container in a logfile ``csobsinfo.log``.
For a direct output on the screen, pass the hidden parameter ``debug=yes`` on the command line:

.. code-block:: bash

  $ csobsinfo debug=yes
  Event list, counts cube, or observation definition file [obs.xml] selected_obs.xml 
  ...
  2016-02-19T15:43:08: +=========+
  2016-02-19T15:43:08: | Summary |
  2016-02-19T15:43:08: +=========+
  2016-02-19T15:43:08:  Unbinned observations .....: 4
  2016-02-19T15:43:08:  Number of events ..........: 6085
  2016-02-19T15:43:08:  Binned observations .......: 0
  2016-02-19T15:43:08:  Number of bins ............: 0
  2016-02-19T15:43:08: 
  2016-02-19T15:43:08: === Pointings ===
  2016-02-19T15:43:08:  Mean zenith angle .........: 46.89
  2016-02-19T15:43:08:  Mean azimuth angle ........: 13.63
  2016-02-19T15:43:08: 
  2016-02-19T15:43:08: === Energy range ===
  2016-02-19T15:43:08:  Emin ......................: 461.48 GeV
  2016-02-19T15:43:08:  Emax ......................: 37.8832 TeV
  2016-02-19T15:43:08: 
  2016-02-19T15:43:08: === Time range ===
  2016-02-19T15:43:08:  Start [MJD] ...............: 53343.9212153
  2016-02-19T15:43:08:  Stop [MJD] ................: 53347.9321065
  2016-02-19T15:43:08:  Total ontime ..............: 6877.00 s = 114.62 min = 1.91 h
  2016-02-19T15:43:08:  Total livetime ............: 6427.52 s = 107.13 min = 1.79 h
  ... 

Further options are described on the reference page :ref:`csobsinfo`. 

Inspecting model XML files
--------------------------
Similar to the observation containers, model XML file can be inspected, too.
The tool :ref:`csmodelinfo` gives a summary of the model container. In particular,
the number of free parameters, number of sky and background models and also number of
parameters that are at their limits might be of interest.

.. code-block:: bash

  $ csmodelinfo debug=yes
  Input model XML file [$CTOOLS/share/models/crab.xml] 
  2016-02-19T15:50:40: +=========+
  2016-02-19T15:50:40: | Summary |
  2016-02-19T15:50:40: +=========+
  2016-02-19T15:50:40: === Instrument specific models ===
  2016-02-19T15:50:40:  All .......................: 1
  2016-02-19T15:50:40:  CTA .......................: 1
  2016-02-19T15:50:40: === Model types ===
  2016-02-19T15:50:40:  PointSource ...............: 1
  2016-02-19T15:50:40:  CTAIrfBackground ..........: 1
  2016-02-19T15:50:40: 
  2016-02-19T15:50:40: +=======================+
  2016-02-19T15:50:40: | Parameter information |
  2016-02-19T15:50:40: +=======================+
  2016-02-19T15:50:40:  All parameters ............: 10
  2016-02-19T15:50:40:  Fixed parameters ..........: 6
  2016-02-19T15:50:40:  Free parameters (total) ...: 4
  2016-02-19T15:50:40:  Free background parameters : 2
  2016-02-19T15:50:40:  Free source parameters ....: 2
  2016-02-19T15:50:40:  Free spectral parameters ..: 2
  2016-02-19T15:50:40:  Parameters at limit .......: 0


Compute upper limit
-------------------
Very often in gamma-ray astronomy sources are at the verge of detection or even not detectable.
In such cases, it is useful to derive an upper limit using :ref:`ctulimit`.

.. code-block:: bash

	$ ctulimit
	Input event list, counts cube or observation definition XML file [selected_obs.xml]
	Input model XML file [crab_models.xml]  
	Source of interest [Crab] 
	
The upper limit will be stored in the log file. To get the limit printed on screen, use the hidden parameter
``debug=yes``.

Compute asymmetric errors
-------------------------
When an analysis approaches its final state, it makes sense to have asymmetric errors on the parameters of the
source of interest. For this purpose, the tool :ref:`cterror` can be used:

.. code-block:: bash

	$ cterror debug=yes
	Input event list, counts cube or observation definition XML file [selected_obs.xml] 
	Input model XML file [crab_models.xml] 
	Source of interest [Crab] 


Compute spectral points
-----------------------

Compute light curves
--------------------

Compute a residual map
----------------------

Compute a test statistics (TS) map
----------------------------------

Compute flux maps
-----------------