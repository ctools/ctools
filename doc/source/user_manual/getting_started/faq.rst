.. _faq:

Frequently asked questions
--------------------------

Below you will find a list of frequently asked questions that has been 
compiled from your user experience and feedback.
In case that you do not find an answer to your question on this page you 
may post your question on the ctools@irap.omp.eu mailing list.
To subscribe to this list you simply need to send an e-mail to
ctools-subscribe@irap.omp.eu (content is irrelevant).
This page will be regularily updated based on the most frequent questions 
asked.

- :ref:`What is the difference between ctools and GammaLib? <faq_ctools_gammalib>`
- :ref:`Are ctools and GammaLib specific to CTA? <faq_instruments>`
- :ref:`I want to upgrade ctools, do I have to upgrade GammaLib? <faq_upgrade>`


.. _faq_ctools_gammalib:

.. topic:: What is the difference between ctools and GammaLib?

   To use a metaphor, ctools is like a village composed of indivdual 
   buildings with different purposes (a school, a church, a townhall, 
   a shop, etc.) while GammaLib is like the bricks that are common to all 
   buildings of the village.

   GammaLib is all about the objects you need to assemble a tool.
   Examples of such objects are "photons", "events", "energies", "times",
   "time intervals", "observations", "pointing directions", etc.
   ctools will put these objects in order to create the functionalities
   that are needed for the analysis of CTA data.
   Examples of such functionalities are "selecting events", "generating a sky 
   map", "fitting a spectrum", "generating a light curve", etc.
   All the inner workings of the objects happen in GammaLib, while ctools 
   essentially provides the concrete to tie the objects together.
   In terms of complexity, 95% of the arithemtics are done in GammaLib while
   ctools is mainly doing some housekeeping.


.. _faq_instruments:

.. topic:: Are ctools and GammaLib specific to CTA?

  The short answer is: ctools yes, GammaLib no.

  ctools have been specifically designed to provide the functionalities that 
  are required to analyse CTA data.
  This does not mean that ctools can not be used for analysing data from 
  other Imaging Air Cherenkov Telescopes; the sole requirement is that the
  data and response files are provided in a CTA compliant format.
  However, ctools are not designed to analyse data from Fermi/LAT, 
  INTEGRAL, COMPTEL or other gamma-ray telescopes.
  Nevertheless, some joint multi-telescope analysis is supported by ctools.

  GammaLib, however, is a universal framework with a core that is completely
  independent of any telescope.
  Plugins are used to interface GammaLib with the data and response files of
  specific telescopes.
  For every new instrument that shall be supported a new plugin needs to be
  developed.
  So far, plugins exists for CTA (and Imaging Air Cherenkov Telescopes in
  general), Fermi/LAT, and COMPTEL.
  Some support for integration of any arbitrary multi-wavelength data 
  exists.


.. _faq_upgrade:

.. topic:: I want to upgrade ctools, do I have to upgrade GammaLib?

  For the moment yes.

  The GammaLib interface is still not fully stabilised, hence ctools is
  continuously adapted to match the interface evolution
  (see the :ref:`download` section for a correspondance of versions).
  We plan however to put the GammaLib interface under change control.
  ctools should then become more independent from GammaLib.
