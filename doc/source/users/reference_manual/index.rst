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
   cslightcrv --- Computes lightcurve <cslightcrv>
   csmodelinfo --- Shows model container content <csmodelinfo>
   csmodelmerge --- Merges several model containers into one file <csmodelmerge>
   csobsdef --- Generates observation definition file <csobsdef>
   csobsinfo --- Shows observation container content <csobsinfo>
   cspull --- Generates pull distribution <cspull>
   csresmap --- Generates residual map <csresmap>
   cssens --- Computes CTA sensitivity <cssens>
   csspec --- Computes spectral points <csspec>
   cstsdist --- Generates TS distribution <cstsdist>
   cstssplit --- Creates commands to slice the ts map computations <cstssplit>
   cstsmapmerge --- Merges slices from ts map computations <cstsmapmerge>
   
 

Scripts to manage an IACT database
----------------------------------

.. toctree::
   :maxdepth: 1

   csobs2caldb --- Creates a caldb entry from an input observation <csobs2caldb>
   csiactdata --- Shows information about IACT data available on the user machine <csiactdata>
   csiactobs --- Generates observation definition file for IACT data from observation IDs <csiactobs>
   csfindobs --- Generates a list of IACT observation IDs <csfindobs>
   csiactcopy --- Copies IACT data from one location to another <csiactcopy>
