.. _start_python:

Using ctools from Python
------------------------

  .. admonition:: What you will learn

     You will learn how to **use the ctools and cscripts from Python** instead
     of typing the commands in the console.

ctools provides two Python modules that allow using all tools and scripts as
Python classes. To use ctools from Python you have to import the ``ctools``
and ``cscripts`` modules into Python. You should also import the ``gammalib``
module, as ctools without GammaLib is generally not very useful.

  .. warning::

     Always import ``gammalib`` before you import ``ctools`` and ``cscripts``.

To simulate the event data from Python you have to type the following (recall
:ref:`how to simulate event data from the command line <start_simulating>` for
comparison) :

.. code-block:: python

   >>> import gammalib
   >>> import ctools
   >>> import cscripts
   >>> sim = ctools.ctobssim()
   >>> sim['inmodel'] = '${CTOOLS}/share/models/crab.xml'
   >>> sim['outevents'] = 'events.fits'
   >>> sim['caldb'] = 'prod2'
   >>> sim['irf'] = 'South_0.5h'
   >>> sim['ra'] = 83.63
   >>> sim['dec'] = 22.01
   >>> sim['rad'] = 5.0
   >>> sim['tmin'] = '2020-01-01T00:00:00'
   >>> sim['tmax'] = '2020-01-01T01:00:00'
   >>> sim['emin'] = 0.03
   >>> sim['emax'] = 200.0
   >>> sim.execute()

The first line after the ``import`` statements generates an instance of the
:ref:`ctobssim` tool as a Python class. User parameters are then set using the
``[ ]`` operator.  After setting all parameters the ``execute()`` method is called to
execute the :ref:`ctobssim` tool. On output the ``events.fits`` FITS file is
created. Until now everything is identical to running the tool from the command
line.

  .. note::

     See the :ref:`reference` for a list of parameters and their types for
     all tools and scripts. You can also inspect all parameters of a specific
     tool by calling the tool follwed by ``--help`` from the command line, e.g.

     .. code-block:: bash

        $ ctobssim --help

Instead of executing the tool you may also *run* the tool using

.. code-block:: python

   >>> sim.run()

The main difference to the ``execute()`` method is that the ``run()`` method
will not write the simulated event file to disk. Why is this useful? Well,
after having typed ``sim.run()`` the :ref:`ctobssim` class still exists as an
object in memory, including all the simulated events. The :ref:`ctobssim`
class has an ``obs()`` method that returns an observation container that holds
the simulated CTA observation with its associated events. To visualise this
container, type:

.. code-block:: python

   >>> print(sim.obs())
   === GObservations ===
    Number of observations ....: 1
    Number of models ..........: 2
    Number of observed events .: 202766
    Number of predicted events : 0

There is one CTA observation in the container and to visualise that observation
type:

.. code-block:: python

   >>> print(sim.obs()[0])
  === GCTAObservation ===
   Name ......................:
   Identifier ................:
   Instrument ................: CTA
   Event file ................: events.fits
   Event type ................: EventList
   Statistics ................: Poisson
   Ontime ....................: 3599.99999976158 s
   Livetime ..................: 3527.99999976635 s
   Deadtime correction .......: 0.98
   User energy range .........: undefined
  === GCTAPointing ===
   Pointing direction ........: (RA,Dec)=(83.63,22.01)
  === GCTAResponseIrf ===
   Caldb mission .............: cta
   Caldb instrument ..........: prod2
   Response name .............: South_0.5h
   Energy dispersion .........: Not used
   Save energy range .........: undefined
  === GCTAEventList ===
   Number of events ..........: 202766 (disposed in "events.fits")
   Time interval .............: 58849.0008007407 - 58849.0424674074 days
   Energy interval ...........: 0.03 - 200 TeV
   Region of interest ........: RA=83.63, DEC=22.01 [0,0] Radius=5 deg

The observation contains a CTA event list that is implement by the GammaLib
class ``GCTAEventList``. You can access the event list using the ``events()``
method. To visualise the individual events you can iterate over the events
using a for loop. This will show the simulated celestial coordinates (RA, DEC),
the coordinate in the camera system [DETX, DETY], the energies and the
terrestrial times (TT) of all events. To do this, type:

.. code-block:: python

   >>> events = sim.obs()[0].events()
   >>> for event in events:
   ...     print(event)
   ...
   Dir=RA=83.6182556152344, DEC=22.2843074798584 [0.00478756989312195,-0.000189682185722814] Energy=63.157357275486 GeV Time=315532804.699509 s (TT)
   Dir=RA=83.8212127685547, DEC=21.9629154205322 [-0.000819836781619686,0.0030950360932479] Energy=53.25997620821 GeV Time=315532805.889628 s (TT)
   Dir=RA=83.6474838256836, DEC=22.0208301544189 [0.000189049359871014,0.000282857291220907] Energy=69.6808248758316 GeV Time=315532806.050809 s (TT)
   ...

You can now benefit from the fact that you have some simulated events in
memory to fit a model to these events using the :ref:`ctlike` class.
You will do this in unbinned mode. Here is what you have to type:

.. code-block:: python

   >>> like = ctools.ctlike(sim.obs())
   >>> like.run()

This is pretty compact. Where are the user parameters?
:ref:`ctlike` doesn't in fact need any parameters as all the relevant
information is already contained in the observation container produced by the
:ref:`ctobssim` class.
And you may have recognised that you constructed the :ref:`ctlike`
instance by using the :ref:`ctobssim` observation container as constructor
argument.

  .. note::

     An observation container, implemented by the ``GObservations`` class
     of GammaLib, is the fundamental brick of any ctools analysis. Many tools
     and scripts handle observation containers, and accept them upon
     construction and return them after running the tool via the ``obs()``
     method.

     Passing observation containers between ctools classes is a very
     convenient and powerful way of building in-memory analysis pipelines.
     However, this implies that you need some computing ressources when
     dealing with large observation containers (for example if you want to
     analyse a few 100 hours of data at once). Also, if the script crashes
     the information is lost.

To check how the fit went you can inspect the optimiser used by
:ref:`ctlike` by typing:

.. code-block:: python

   >>> print(like.opt())         
   === GOptimizerLM ===
    Optimized function value ..: 757618.425
    Absolute precision ........: 0.005
    Acceptable value decrease .: 2
    Optimization status .......: converged
    Number of parameters ......: 10
    Number of free parameters .: 4
    Number of iterations ......: 2
    Lambda ....................: 1e-05

You see that the fit converged after 2 iterations. Out of 10 parameters in the
model 4 have been fitted (the others were kept fixed). To inspect the fit
results you can print the model container that can be access using the
``models()`` method of the observation container:

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
    Model components ..........: "PointSource" * "PowerLaw" * "Constant"
    Number of parameters ......: 6
    Number of spatial par's ...: 2
     RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
     DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
    Number of spectral par's ..: 3
     Prefactor ................: 5.73901839854211e-16 +/- 5.75150678022232e-18 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
     Index ....................: -2.46442704108736 +/- 0.00738985717701118 [-0,-5]  (free,scale=-1,gradient)
     PivotEnergy ..............: 300000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
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
     Prefactor ................: 0.995734795709415 +/- 0.00725726198557001 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
     Index ....................: -0.00198344390952578 +/- 0.00257263793207067 [-5,5]  (free,scale=1,gradient)
     PivotEnergy ..............: 1000000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
    Number of temporal par's ..: 1
     Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

Suppose you want to repeat the fit by optimising also the position of the 
point source. This is easy from Python. Type the following:

.. code-block:: python

   >>> like.obs().models()['Crab']['RA'].free()
   >>> like.obs().models()['Crab']['DEC'].free()
   >>> like.run()
   >>> print(like.obs().models())
   ...
     RA .......................: 83.6334964430607 +/- 0.000614768905318078 [-360,360] deg (free,scale=1)
     DEC ......................: 22.0149467750683 +/- 0.000559548446961086 [-90,90] deg (free,scale=1)


The ``like.obs().models()`` method provides the model container, using the 
``['Crab']`` operator you access the Crab model in that container and using
the ``['RA']`` and ``['DEC']`` operators you access the relevant model
parameters. The ``free()`` method frees a parameter, the opposite would be a
call to the ``fix()`` method.
