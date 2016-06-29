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
- :ref:`Why does my ctool complain about missing parameters? <faq_pars>`
- :ref:`Should I used binned or unbinned analysis? <faq_analysis>`
- :ref:`How precise are ctools? <faq_precision>`


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


.. _faq_pars:

.. topic:: Why does my ctool complain about missing parameters?

  After upgrading to a new ctools version, a ctool may issue the following
  error

  ``*** ERROR in GApplicationPars::operator[](std::string&): Invalid 
  argument. Parameter "XXX" has not been found in parameter file.
  Please specify a valid parameter name.``

  This error may occur if some parameters have been added to a ctool and
  you still have the old parameter file without that parameter sitting in
  your pfiles folder. Simply erasing all files in the pfiles folder should
  fix this problem.


.. _faq_analysis:

.. topic:: Should I used binned or unbinned analysis?

  This depends on the amount of data you want to analyse and to some extent
  on the question you want to answer with your analysis.
  For an unbinned analysis, the computation time increases about linearly
  with the number of events and hence with the duration of the observation,
  while for binned analysis the computation time depends only on the number
  of bins that is used.
  Consequently, for short observation times (below 30 hours), unbinned
  analysis is faster, while for longer times it is advantageous to use
  a binned analysis.
  If you go to very low energies, either make sure that you have enough
  energy bins in a binned analysis to sample properly the strong drop in the
  effective area towards small energies, or better, use unbinned analysis.
  If you'd like to fit a Gaussian line to your data you may also prefer
  unbinned over binned analysis, as the fine sampling required to resolve
  the line may require a prohibitive large number of energy bins for a binned
  analysis.


.. _faq_precision:

.. topic:: How precise are ctools?

  The ctools and gammalib codes have a numerical accuracy of better than 1%.
  This means that if you use ctools to determine for example the flux 
  received from a source or the spectral points of an SED, the relative 
  precision of the flux or the spectral points is better than 1%.
  The same is true for spatial parameters, such as source position or
  source extension.
  For many cases the actual numerical precision is in fact much better
  than 1%, but in any case, it should never be worse.
  Note, however, that this does not imply that source parameters can be
  determined with CTA with an accuracy of 1%. The accuracy depends in the
  end on the precision to which the instrument response function is known,
  which should be more in the 10% range.
