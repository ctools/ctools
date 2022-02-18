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

