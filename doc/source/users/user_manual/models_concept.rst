.. _sec_models_concept:

Overall concept
------------------

A central element of any data analysis is the modeling of the observed
event distribution.
Generally, the measured events can be separated into two distinct classes:
events being attributed to gamma rays of celestial origin (the "source" 
events) or events being attributed to any other kind of trigger (the 
"background" events).
In the first case, the :ref:`Instrument Response Functions (IRFs) <um_response>`
describe how an incident gamma-ray converts into a measured event.
In the second case there is often no general prescription, and the 
distribution of background events is commonly modeled directly in the data 
space of the observable quantities.
The main difficulty of gamma-ray astronomy is that not all events can be
tagged a priori as being source or background, hence their separation can 
only be done on a statistical basis, based on the different morphologies, 
spectral characteristics or eventually temporal signatures of both event
categories.

For this purpose, ctools use a general model that describes the spatial, 
spectral and temporal properties of the source and background components.
The model is composed of an arbitrary number of model components that
add up linearly to provide a prediction of the expected number of events.
Model components are generally parametric, and model parameters can be 
adjusted through a maximum likelihood procedure to find the set of 
parameters that represent best the measured data.
Model components representing celestial sources are convolved with the 
:ref:`IRFs <um_response>` to predict the expected number of source events in
the data.
Background model components will be directly expressed as expected number 
of background events without any :ref:`IRF <um_response>` convolution.
