.. _reference:

Reference Manual
================

This manual provides reference information for all ctools and csripts.
General information on ctools usage can be found `here
<../user_manual/introduction.html>`_. The description of user parameters in the reference
manual is documented `here <usage.html>`__.


Below you find links to the command line reference for all available tools and scripts.

Generic analysis tools and scripts
----------------------------------

.. toctree::
   :maxdepth: 1

   ctbutterfly --- Compute butterfly <ctbutterfly>
   cterror --- Calculates likelihood profile errors <cterror>
   ctlike --- Performs maximum likelihood fitting <ctlike>
   ctmapcube --- Generates a map cube <ctmapcube>
   cttsmap --- Generates Test Statistic map <cttsmap>
   ctulimit --- Calculates upper limit <ctulimit>
   cscaldb --- Lists available instrument response functions <cscaldb>
   cslightcrv --- Computes light curve <cslightcrv>
   csmodelinfo --- Shows model container content <csmodelinfo>
   csmodelmerge --- Merges several model containers into one file <csmodelmerge>
   csmodelselect --- Select models from model definition file <csmodelselect>
   csmodelsois --- Generate map cube from subset of models <csmodelsois>
   csobsinfo --- Shows observation container content <csobsinfo>
   csspec --- Computes spectral points <csspec>
   cstsdist --- Generates Test Statistic distribution <cstsdist>
   cstsmapmerge --- Merges slices from Test Statistic map computations <cstsmapmerge>
   cstsmapsplit --- Creates commands to split the Test Statistic map computations <cstsmapsplit>


CTA and IACT analysis tools and scripts
---------------------------------------

.. toctree::
   :maxdepth: 1

   ctbin --- Generates counts cube <ctbin>
   ctcubemask --- Filter counts cube <ctcubemask>
   ctexpcube --- Generates exposure cube <ctexpcube>
   ctpsfcube --- Generates point spread function cube <ctpsfcube>
   ctedispcube --- Generates energy dispersion cube <ctedispcube>
   ctbkgcube --- Generates background cube <ctbkgcube>
   ctfindvar --- Search for source variability <ctfindvar>
   ctmodel --- Computes model counts cube <ctmodel>
   ctobssim --- Simulate observations <ctobssim>
   ctphase --- Computes the phase of each event <ctphase>
   ctprob --- Computes event probability for a given model <ctprob>
   ctselect --- Selects event data <ctselect>
   ctskymap --- Generates sky map <ctskymap>
   csbkgmodel --- Generates background model for 3D analysis <csbkgmodel>
   csebins --- Generates energy boundaries for stacked analysis <csebins>
   csobsdef --- Generates observation definition file <csobsdef>
   csobsselect --- Select observations from observation definition file <csobsselect>
   csphagen --- Generates PHA, ARF, RMF files based on source/background regions <csphagen>
   csphasecrv --- Computes phase curve <csphasecrv>
   cspull --- Generates pull distribution <cspull>
   csresmap --- Generates residual map <csresmap>
   csresspec --- Generates residual spectrum <csresspec>
   csadd2caldb --- Adds CTA response function to calibration database <csadd2caldb>
   csroot2caldb --- Creates a caldb entry from a ROOT file <csroot2caldb>
   csscs --- Performs spectral component separation <csscs>
   cssens --- Computes CTA sensitivity <cssens>
   cssrcdetect --- Detects sources in sky map <cssrcdetect>
   csviscube --- Computes visibility cube <csviscube>


IACT database management scripts
--------------------------------

.. warning::

   The ``csiactdata``, ``csiactobs``, ``csfindobs`` and ``csiactcopy`` scripts
   rely on the ``json`` Python module which is only available in Python 2.6
   or higher. These scripts will not work on older Python versions.

.. toctree::
   :maxdepth: 1

   csobs2caldb --- Creates a caldb entry from an input observation <csobs2caldb>
   csiactdata --- Shows information about IACT data available on the user machine <csiactdata>
   csiactobs --- Generates observation definition file for IACT data from observation IDs <csiactobs>
   csfindobs --- Generates a list of IACT observation IDs <csfindobs>
   csiactcopy --- Copies IACT data from one location to another <csiactcopy>
   

COMPTEL Science Analysis scripts
--------------------------------

.. toctree::
   :maxdepth: 1

   comlixfit --- Fit model to data using SRCLIX algorithm <comlixfit>
   comlixmap --- Create TS map using SRCLIX algorithm <comlixmap>
   comobsadd --- Combine observations <comobsadd>
   comobsback --- Generate background model for COMPTEL observations <comobsback>
   comobsbin --- Bin COMPTEL observations <comobsbin>
   comobsmodel --- Generate model for binned COMPTEL observations <comobsmodel>
   comobsres --- Generate residuals of COMPTEL observations <comobsres>
   comobsselect --- Select observations from COMPTEL database <comobsselect>
   comobssim --- Simulate COMPTEL observations <comobssim>
   compulbin --- Generate pulse profiles for pulsars <compulbin>
   comsrcdetect --- Detect source in TS map <comsrcdetect>
