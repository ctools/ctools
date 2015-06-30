.. _python:

Using ctools from Python
------------------------

Getting started
~~~~~~~~~~~~~~~

ctools provides a Python module that allows using all tools and scripts as 
Python classes.
To use ctools from Python all you have to do is to import the ctools 
module into Python.
You should also import the GammaLib module, as ctools without GammaLib is 
generally not very useful.

.. code-block:: python

   >>> import gammalib
   >>> import ctools

To illustrate how to use ctools from Python, below is a working example 
of an event list simulation using the :ref:`ctobssim` class.
An instance ``sim`` of the :ref:`ctobssim` class is generated and user 
parameters are set for this instance using the ``[ ]`` operator.
Note that specific methods exist for [file], [string], [real], [integer], 
and [boolean] parameter types (see the :ref:`reference` for a list of 
parameters and their types).
The ``execute()`` method executes the :ref:`ctobssim` class in the same 
way as if it were executed from the command line.

.. code-block:: python

   >>> import gammalib
   >>> import ctools
   >>> sim = ctools.ctobssim()
   >>> sim["inmodel"].filename("${CTOOLS}/share/models/crab.xml")
   >>> sim["outevents"].filename("events.fits")
   >>> sim["caldb"].string("prod2")
   >>> sim["irf"].string("South_50h")
   >>> sim["ra"].real(83.63)
   >>> sim["dec"].real(22.01)
   >>> sim["rad"].real(5.0)
   >>> sim["tmin"].real(0.0)
   >>> sim["tmax"].real(1800.0)
   >>> sim["emin"].real(0.1)
   >>> sim["emax"].real(100.0)
   >>> sim.execute()

Alternatively, you may "run" the :ref:`ctobssim` tool using

.. code-block:: python

   >>> sim.run()

The main difference to the ``execute()`` method is that the ``run()`` 
will not write the simulated event file to disk.
Why is this useful?
Well, after having typed ``sim.run()`` the :ref:`ctobssim` class still 
exists as an object in memory, including all the simulated events.
The :ref:`ctobssim` class has an ``obs()`` method that returns an 
observation container that holds the simulated CTA observation with its 
associated events.
To visualise this container, type:

.. code-block:: python

   >>> print(sim.obs())
   === GObservations ===
    Number of observations ....: 1
    Number of predicted events : 0
   === GCTAObservation ===
    Name ......................: 
    Identifier ................: 
    Instrument ................: CTA
    Event file ................: events.fits
    Event type ................: EventList
    Statistics ................: Poisson
    Ontime ....................: 1800 s
    Livetime ..................: 1710 s
    Deadtime correction .......: 0.95
    User energy range .........: undefined
   ...

There is one CTA observation in the container and to visualise the events 
in that observation you may type:

.. code-block:: python

   >>> print(sim.obs()[0].events())
   === GCTAEventList ===
    Number of events ..........: 31928
    Time interval .............: 51544.5 - 51544.5 days
   === GEbounds ===
    Number of intervals .......: 1
    Energy range ..............: 100 GeV - 100 TeV
   === GCTARoi ===
    ROI centre ................: RA=83.63, DEC=22.01 [0,0]
    ROI radius ................: 5 deg

The ``obs()[0]`` operator returns the first observation in the observation 
container, the ``events()`` operator returns the event list in that 
observation.
To see what kind of object you actually got, use:

.. code-block:: python

   >>> type(sim.obs()[0].events())
   <class 'gammalib.cta.GCTAEventList'>

The CTA event list is implement by the ``GCTAEventList`` class in the 
``cta`` module of GammaLib.
To visualise the individual events you can iterate over the events using a 
for loop.
This will show the simulated celestial coordinates (RA, DEC), the 
coordinate in the camera system [DETX, DETY], energies and 
terrestrial times (TT) of all events. 

.. code-block:: python

   >>> events = sim.obs()[0].events()
   >>> for event in events:
   ...     print(event)
   ...
   Dir=RA=83.5759, DEC=21.9757 [-0.000597839,-0.000875483] Energy=328.153 GeV Time=-3.15576e+08 s (TT)
   Dir=RA=83.5775, DEC=22.0128 [4.88943e-05,-0.000849504] Energy=353.264 GeV Time=-3.15576e+08 s (TT)
   Dir=RA=83.4302, DEC=21.931 [-0.00137616,-0.0032343] Energy=117.999 GeV Time=-3.15576e+08 s (TT)
   ...

.. note::

   If you inspect the event list in detail you may recognise that all 
   events appear twice in the event list.
   This is because we ran :ref:`ctobssim` twice (a first time when calling
   the ``execute()`` method and a second time when calling the ``run()`` 
   method).
   **ctobssim will append simulated events to any pre-existing
   events in an observation container.**
   To get a single simulation you should clear the observation container 
   and re-run :ref:`ctobssim`:

   .. code-block:: python

     >>> sim.obs().clear()           
     >>> sim.run()                   
     >>> print(sim.obs()[0].events())
     === GCTAEventList ===
      Number of events ..........: 15964
      Time interval .............: 51544.5 - 51544.5 days
     === GEbounds ===
      Number of intervals .......: 1
      Energy range ..............: 100 GeV - 100 TeV
     === GCTARoi ===
      ROI centre ................: RA=83.63, DEC=22.01 [0,0]
      ROI radius ................: 5 deg


Performing a maximum likelihood analysis in Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can now benefit from the fact that we have some simulated events in 
memory to fit a model to these events using the :ref:`ctlike` class.
We will do this in unbinned mode.
Here is what you have to do:

.. code-block:: python

   >>> like = ctools.ctlike(sim.obs())
   >>> like.run()

This is pretty compact.
Where are the user parameters?
:ref:`ctlike` doesn't in fact need any as all the relevant information is 
already contained in the observation container produced by the 
:ref:`ctobssim` class.
And you may have recognised that we constructed the :ref:`ctlike` 
instance by using the :ref:`ctobssim` observation container as
constructor argument.

To check how the fit went you can inspect the optimiser used by
:ref:`ctlike`:

.. code-block:: python

   >>> print(like.opt())         
   === GOptimizerLM ===
    Optimized function value ..: 101048.382
    Absolute precision ........: 0.005
    Acceptable value decrease .: 2
    Optimization status .......: converged
    Number of parameters ......: 10
    Number of free parameters .: 4
    Number of iterations ......: 2
    Lambda ....................: 1e-05

Apparently, the fit converged after 2 iterations.
Out of 10 parameters in the model 4 have been fitted (the others were kept
fixed).
To inspect the fit results you can print the model container that is a 
member of the observation container:

.. code-block:: python

   >>> print(like.obs().models())
   === GModels ===
    Number of models ..........: 2
    Number of parameters ......: 10
   === GModelSky ===
    Name ......................: Crab
    Instruments ...............: all
    Instrument scale factors ..: unity
    Observation identifiers ...: all
    Model type ................: PointSource
    Model components ..........: "SkyDirFunction" * "PowerLaw" * "Constant"
    Number of parameters ......: 6
    Number of spatial par's ...: 2
     RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
     DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
    Number of spectral par's ..: 3
     Prefactor ................: 5.80565e-16 +/- 1.11147e-17 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
     Index ....................: -2.49266 +/- 0.0163476 [-0,-5]  (free,scale=-1,gradient)
     PivotEnergy ..............: 300000 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
    Number of temporal par's ..: 1
     Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   === GCTAModelIrfBackground ===
    Name ......................: CTABackgroundModel
    Instruments ...............: CTA
    Instrument scale factors ..: unity
    Observation identifiers ...: all
    Model type ................: "PowerLaw" * "Constant"
    Number of parameters ......: 4
    Number of spectral par's ..: 3
     Prefactor ................: 1.01701 +/- 0.0216563 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
     Index ....................: 0.0230864 +/- 0.0123006 [-5,5]  (free,scale=1,gradient)
     PivotEnergy ..............: 1e+06 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
    Number of temporal par's ..: 1
     Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

Suppose you want to repeat the fit by optimising also the position of the 
point source.
This is easy from Python:

.. code-block:: python

   >>> like.obs().models()["Crab"]["RA"].free()
   >>> like.obs().models()["Crab"]["DEC"].free()
   >>> like.run()
   >>> print(like.obs().models())
   ...
     RA .......................: 83.6341 +/- 0.000958947 [-360,360] deg (free,scale=1)
     DEC ......................: 22.0139 +/- 0.000870912 [-90,90] deg (free,scale=1)

The ``like.obs().models()`` method provides the model container, using the 
``["Crab"]`` operator we access the Crab model in that container and using 
the ``["RA"]`` and ``["DEC"]`` methods we access the relevant model 
parameters.
The ``free()`` method frees a parameter, the opposite would be a call to 
the ``fix()`` method.

.. warning::

   Passing observation containers between ctools classes is a very 
   convenient and powerful way of building in-memory analysis pipelines.
   However, this implies that you need some computing ressources when 
   dealing with large observation containers (for example if you want to 
   analyse a few 100 hours of data at once).

.. warning::

   You have to be aware about the scope of the objects you're 
   dealing with.
   In the above example, the ``sim.obs()`` container is allocated by the
   :ref:`ctobssim` class, hence it disappears (a.k.a. goes out of scope)
   once the :ref:`ctobssim` class is deleted, as illustrated by the 
   following example:

   .. code-block:: python

      >>> obs = sim.obs()
      >>> del sim
      >>> print(obs)
      Segmentation fault

   Note that

   .. code-block:: python

      >>> obs = sim.obs()
      >>> del sim
      >>> print(like.obs())
      === GObservations ===
       Number of observations ....: 1
       Number of predicted events : 15964

   is okay as the ``ctools.ctlike(sim.obs())`` constructor will create
   a copy of the observation container that lives within the :ref:`ctlike`
   instance.
   To preserve an observation container after a ctools object goes out 
   of scope you have to create a local copy of the container using the
   ``copy()`` method:

   .. code-block:: python

      >>> obs = sim.obs().copy()
      >>> del sim
      >>> print(obs)
      === GObservations ===
       Number of observations ....: 1
       Number of predicted events : 0


Using obsutils
~~~~~~~~~~~~~~

ctools provides the Python module ``obsutils`` that may further simplify 
your analysis efforts.
``obsutils`` is a Python script that makes use of the GammaLib and ctools 
modules to create standard analysis steps.
Here an example of how to use ``obsutils``:

.. code-block:: python

   >>> import gammalib
   >>> import ctools
   >>> from ctools import obsutils
   >>> pattern = obsutils.set_obs_patterns("four", ra=83.63, dec=22.01, offset=1.0)
   >>> obs = obsutils.set_obs_list(pattern, duration=1800, emin=0.1, emax=100.0, rad=5.0, caldb="prod2", irf="South_50h")
   >>> print(obs)   
   === GObservations ===
    Number of observations ....: 4
    Number of predicted events : 0
   ...
   >>> obs.models(gammalib.GModels("${CTOOLS}/share/models/crab.xml"))
   >>> obs = obsutils.sim(obs)
   >>> like = obsutils.fit(obs)
   >>> print(like.obs().models())   
   === GModels ===
    Number of models ..........: 2
    Number of parameters ......: 10
   ...

The module is imported using the ``from ctools import obsutils`` directive.
The ``obsutils.set_obs_patterns()`` function will create a pointing 
pattern of four observations located at offset angles of 1 degree from the 
nominal location of the Crab nebula.
The ``obsutils.set_obs_list()`` will build an observation container from 
that pattern where each pointing will have a duration of 1800 seconds, 
cover the 0.1-100 TeV energy range and a field of view of 5 degrees.
The standard dummy CTA calibration information will be used.
A model is then appended to the observation container using the
``obs.models()`` method.
The ``obsutils.sim()`` function then simulates the event data, the
``obsutils.fit()`` function performs a maximum likelihood fit.

.. note::

   The ``obsutils`` module is not yet fully developed and more convenience 
   functions will be added in the future.


Beyond the first steps
~~~~~~~~~~~~~~~~~~~~~~

You now have learned the basics of using ctools and GammaLib within Python.
To go beyond these initial steps you may check the Python scripts in the
``examples`` folder that provide useful analysis examples.
Check the ``README`` file in that folder for an explanation of the scripts.


