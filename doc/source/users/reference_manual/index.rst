.. _reference:

Reference Manual
================

This manual provides reference information for all ctools and csripts.
General information on ctools usage can be found `here <usage.html>`__.

Below you find links to the command line reference for the tools and scripts
that are available.

ctools
------

.. toctree::
   :maxdepth: 1

   ctbin --- Generates counts cube <ctbin>
   ctbkgcube --- Generates background cube <ctbkgcube>
   ctbutterfly --- Compute butterfly <ctbutterfly>
   ctcubemask --- Filter counts cube <ctcubemask>
   ctedispcube --- Generates energy dispersion cube <ctedispcube>
   cterror --- Calculates likelihood profile errors <cterror>
   ctexpcube --- Generates exposure cube <ctexpcube>
   ctlike --- Performs maximum likelihood fitting <ctlike>
   ctmapcube --- Generates a map cube <ctmapcube>
   ctmodel --- Computes model counts cube <ctmodel>
   ctobssim --- Simulate observations <ctobssim>
   ctpsfcube --- Generates point spread function cube <ctpsfcube>
   ctselect --- Selects event data <ctselect>
   ctskymap --- Generates sky map <ctskymap>
   cttsmap --- Generates Test Statistic map <cttsmap>
   ctulimit --- Calculates upper limit <ctulimit>


cscripts
--------

.. toctree::
   :maxdepth: 1

   cscaldb --- Lists available instrument response functions <cscaldb>
   csebins --- Generates energy boundaries for stacked analysis <csebins>
   cslightcrv --- Computes lightcurve <cslightcrv>
   csmodelinfo --- Shows model container content <csmodelinfo>
   csmodelmerge --- Merges several model containers into one file <csmodelmerge>
   csmodelselect --- Select models from model definition file <csmodelselect>
   csobsdef --- Generates observation definition file <csobsdef>
   csobsinfo --- Shows observation container content <csobsinfo>
   csobsselect --- Select observations from observation definition file <csobsselect>
   cspull --- Generates pull distribution <cspull>
   csresmap --- Generates residual map <csresmap>
   cssens --- Computes CTA sensitivity <cssens>
   csspec --- Computes spectral points <csspec>
   cssrcdetect --- Detects sources in sky map <cssrcdetect>
   cstsdist --- Generates Test Statistic distribution <cstsdist>
   cstsmapsplit --- Creates commands to split the Test Statistic map computations <cstsmapsplit>
   cstsmapmerge --- Merges slices from Test Statistic map computations <cstsmapmerge>
   csviscube --- Computes visibility cube <csviscube>


Scripts to manage an IACT database
----------------------------------

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
