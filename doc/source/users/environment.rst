.. _sec_environment:

Environmental impact
====================

Like any human activity, using ctools will have an environmental impact. This
includes for example the greenhouse gas emissions related to electricity
production that is needed to run your computer or data centre. It includes as
well the environmental footprint of your computing hardware, the power and
material behind your internet connection, your storage array, and even the
footprint of the personnel that eventually runs your data centre.

As ctools user you should be aware of this footprint, and you should aim
to reduce it as much as possible.

To guide you in this endeavor, all ctools and cscripts inform you at the end
of the log file about the estimated carbon footprint of your computations:

.. code-block:: bash

   2022-02-18T22:33:26: Application "ctlike" terminated after 3910.0 wall clock seconds, consuming 3881.18 seconds of CPU time and generating a carbon footprint of 3.81083 g eCO2.

This estimate is based on the study of
`Berthoud et al. (2020) <https://hal.archives-ouvertes.fr/hal-02549565v4/document>`_
who assessed the carbon footprint of the DAHU cluster of the UMS GRICAD (Grenoble,
France) in 2019. Their study included

* server fabrication
* server environment
* server usage (electricity)
* travelling of personnel in the context of their work
* travelling of personnel from home to office
* personnel equipment
* personnel energy

`Berthoud et al. (2020) <https://hal.archives-ouvertes.fr/hal-02549565v4/document>`_
found that, per CPU hour, the footprint of the GRICAD data centre was 2.43 g CO2e
(equivalent CO2 emissions) due to electricity use and 2.25 g CO2e due to the data
centre environment, hardware fabrication and the personnel. For their study a
carbon intensity of 108 g eCO2 / kWh for greenhouse gas emission due to
electricity generation was assumed.

The latter varies substantially between countries. While the current electricity mix
in France has a (rounded) carbon intensity of 60 g eCO2 / kWh, values go through
240 g eCO2 / kWh (Spain), 410 g eCO2 / kWh (Italy), 460 g eCO2 / kWh (Germany),
520 g eCO2 / kWh (USA), up to 840 g eCO2 / kWh (Australia). So eventually, you should
think about doing your computations in a country with a low carbon intensity, or make
sure that your data centre is running on a low carbon energy source.

ctools will determine in which country your code is executed and fetch the appropriate
carbon intensity for the general energy mix in this country to estimate the carbon
footprint of your computations. The estimate is done using

.. math::
   {\rm Carbon\,\,footprint} = \left( 2.25 + 2.43 \frac{\rm CI}{108} \right) \times {\rm CPU\,\,hours}

where :math:`{\rm CI}` is the carbon intensity of the electricity mix in your country
in units of g eCO2 / kWh, :math:`{\rm CPU\,\,hours}` is the number of CPU hours your
computation took, and :math:`{\rm Carbon\,\,footprint}` is the carbon footprint in g eCO2.

The carbon footprint of all ctools and cscript runs are stored in an ASCII file
that you will find in a folder under your home directory, named ``.gamma/statistics.csv``.
This file is regularly scanned and purged by the GammaLib daemon that creates the
high-level statistics file ``.gamma/statistics.xml`` that is in the same folder. You may
display the content of this file using the ``csfootprint`` script:

.. code-block:: bash

   $ csfootprint debug=yes
   Start time for report (UTC string, JD, MJD or MET in seconds) [NONE]
   Stop time for report (UTC string, JD, MJD or MET in seconds) [NONE]
   Output graphics file (NONE if no graphics should be generated) [NONE]
   ...
   2022-02-20T14:38:29: +===================+
   2022-02-20T14:38:29: | Global statistics |
   2022-02-20T14:38:29: +===================+
   2022-02-20T14:38:29:  Creation date .............: 2022-02-19T22:09:05
   2022-02-20T14:38:29:  Last statistics update ....: 2022-02-20T14:35:20
   2022-02-20T14:38:29:  Statistics date interval ..: 2022-02-18T22:21:56 - 2022-02-20T14:35:17
   2022-02-20T14:38:29:  Used date interval ........: 2022-02-18T22:21:56 - 2022-02-20T14:35:17
   2022-02-20T14:38:29:  Duration of used interval .: 40.222 hours
   2022-02-20T14:38:29:  Total number of ctool runs : 3104
   2022-02-20T14:38:29:  Total wall clock time .....: 34.200 minutes
   2022-02-20T14:38:29:  Total CPU time ............: 27.543 minutes
   2022-02-20T14:38:29:  Average CPU load ..........: 80.5 %
   2022-02-20T14:38:29:  Total carbon footprint ....: 1.623 g CO2e
   2022-02-20T14:38:29:  Average carbon intensity ..: 3.535 g CO2e / CPU hour
   2022-02-20T14:38:29:  Average daily footprint ...: 0.968 g CO2e / day
   2022-02-20T14:38:29:  Expected annual footprint .: 353.632 g CO2e / year
   2022-02-20T14:38:29:
   2022-02-20T14:38:29: +==================+
   2022-02-20T14:38:29: | Daily statistics |
   2022-02-20T14:38:29: +==================+
   2022-02-20T14:38:29: === Carbon footprint ===
   2022-02-20T14:38:29:  2022-02-18 ................: 1.336 g CO2e
   2022-02-20T14:38:29:  2022-02-19 ................: 0.010 g CO2e
   2022-02-20T14:38:29:  2022-02-20 ................: 0.277 g CO2e
   2022-02-20T14:38:29: === ctools or cscript calls ===
   2022-02-20T14:38:29:  2022-02-18 ................: 2895
   2022-02-20T14:38:29:  2022-02-19 ................: 80
   2022-02-20T14:38:29:  2022-02-20 ................: 129
   2022-02-20T14:38:29: === Used wall clock time ===
   2022-02-20T14:38:29:  2022-02-18 ................: 25.483 minutes
   2022-02-20T14:38:29:  2022-02-19 ................: 1.033 minutes
   2022-02-20T14:38:29:  2022-02-20 ................: 7.683 minutes
   2022-02-20T14:38:29: === Used CPU time ===
   2022-02-20T14:38:29:  2022-02-18 ................: 22.682 minutes
   2022-02-20T14:38:29:  2022-02-19 ................: 9.895 seconds
   2022-02-20T14:38:29:  2022-02-20 ................: 4.695 minutes
   2022-02-20T14:38:29:
   2022-02-20T14:38:29: +================================+
   2022-02-20T14:38:29: | ctools and cscripts statistics |
   2022-02-20T14:38:29: +================================+
   2022-02-20T14:38:29: === Carbon footprint ===
   2022-02-20T14:38:29:  ctlike ....................: 0.318 g CO2e
   2022-02-20T14:38:29:  cspull ....................: 0.186 g CO2e
   2022-02-20T14:38:29:  csbkgmodel ................: 0.178 g CO2e
   2022-02-20T14:38:29:  ctobssim ..................: 0.172 g CO2e
   2022-02-20T14:38:29:  ctobssim ..................: 0.087 g CO2e
   2022-02-20T14:38:29:  comobsbin .................: 0.081 g CO2e
   2022-02-20T14:38:29:  comlixmap .................: 0.069 g CO2e
   2022-02-20T14:38:29:  ctulimit ..................: 0.064 g CO2e
   2022-02-20T14:38:29:  cstsdist ..................: 0.058 g CO2e
   2022-02-20T14:38:29:  ctselect ..................: 0.044 g CO2e
   2022-02-20T14:38:29:  ... (list truncated after 10 entries) ...
   2022-02-20T14:38:29: === ctools or cscript calls ===
   2022-02-20T14:38:29:  ctobssim ..................: 1795
   2022-02-20T14:38:29:  ctlike ....................: 220
   2022-02-20T14:38:29:  ctselect ..................: 87
   2022-02-20T14:38:29:  csfootprint ...............: 80
   2022-02-20T14:38:29:  cscaldb ...................: 79
   2022-02-20T14:38:29:  csscs .....................: 55
   2022-02-20T14:38:29:  ctulimit ..................: 54
   2022-02-20T14:38:29:  ctcubemask ................: 51
   2022-02-20T14:38:29:  csphagen ..................: 46
   2022-02-20T14:38:29:  ctbin .....................: 31
   2022-02-20T14:38:29:  ... (list truncated after 10 entries) ...
   2022-02-20T14:38:29: === Used wall clock time ===
   2022-02-20T14:38:29:  ctlike ....................: 5.367 minutes
   2022-02-20T14:38:29:  cspull ....................: 3.817 minutes
   2022-02-20T14:38:29:  csbkgmodel ................: 3.067 minutes
   2022-02-20T14:38:29:  ctobssim ..................: 2.967 minutes
   2022-02-20T14:38:29:  csfootprint ...............: 2.300 minutes
   2022-02-20T14:38:29:  ctobssim ..................: 1.483 minutes
   2022-02-20T14:38:29:  comobsbin .................: 1.367 minutes
   2022-02-20T14:38:29:  cstsdist ..................: 1.300 minutes
   2022-02-20T14:38:29:  comlixmap .................: 1.183 minutes
   2022-02-20T14:38:29:  ctselect ..................: 1.083 minutes
   2022-02-20T14:38:29:  ... (list truncated after 10 entries) ...
   2022-02-20T14:38:29: === Used CPU time ===
   2022-02-20T14:38:29:  ctlike ....................: 5.392 minutes
   2022-02-20T14:38:29:  cspull ....................: 3.157 minutes
   2022-02-20T14:38:29:  csbkgmodel ................: 3.029 minutes
   2022-02-20T14:38:29:  ctobssim ..................: 2.921 minutes
   2022-02-20T14:38:29:  ctobssim ..................: 1.469 minutes
   2022-02-20T14:38:29:  comobsbin .................: 1.367 minutes
   2022-02-20T14:38:29:  comlixmap .................: 1.169 minutes
   2022-02-20T14:38:29:  ctulimit ..................: 1.086 minutes
   2022-02-20T14:38:29:  cstsdist ..................: 58.978 seconds
   2022-02-20T14:38:29:  ctselect ..................: 45.005 seconds
   2022-02-20T14:38:29:  ... (list truncated after 10 entries) ...


