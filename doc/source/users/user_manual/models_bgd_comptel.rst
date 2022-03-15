.. _um_models_bgd_comptel:

COMPTEL background models
-------------------------

The following section presents the model that is available in ctools for the
modelling of instrumental background in COMPTEL data.

In general, the instrumental background in COMPTEL data for a given energy
band :math:`E'` is specified by a so-called DRB cube:

.. math::
   M(p',E') = {\rm DRB}(\chi, \psi, \bar{\varphi} | E')

where :math:`M(p',E')` is given in units of
:math:`{\rm events} \,\, {\rm s}^{-1} {\rm MeV}^{-1} {\rm sr}^{-1}`.
:math:`(\chi, \psi)` is the photon scatter direction while
:math:`\bar{\varphi}` is the Compton scattering angle.


DRB cubes
^^^^^^^^^

DRB cubes are specified as FITS files in the
:ref:`observation definition file <glossary_obsdef>`
using the ``DRB`` parameter, e.g.

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <observation_list title="observation library">
     <observation name="Crab" id="100001" instrument="COM">
       <parameter name="DRE" file="m50438_dre.fits"/>
       <parameter name="DRB" file="bgdlix_drb.fits"/>
       <parameter name="DRG" file="m34997_drg.fits"/>
       <parameter name="DRX" file="m32171_drx.fits"/>
       <parameter name="IAQ" value="SIM2(0.75-1.00)MeV(2)deg"/>
     </observation>
     ...
   </observation_list>

DRB files are not provided by HEASARC. They may be generated from the data
using the :ref:`comobsback` script.


DRB fitting
^^^^^^^^^^^

The usual way of handling the COMPTEL background model during a maximum likelihood
fit is to adjust the normalisation of each :math:`\bar{\varphi}`-layer of the
DRB cube to the data. This is accomplished by using the ``DRBPhibarBins`` model.
Below is an example :ref:`model definition file <glossary_moddef>` for this
model.

.. code-block:: xml

   <source name="Background(0.75-1.0)MeV" type="DRBPhibarBins" instrument="COM" id="100001">
     <parameter name="Normalization" value="1" scale="1" min="0" max="1000" free="0" />
     <parameter name="Normalization" value="1" scale="1" min="0" max="1000" free="0" />
     <parameter name="Normalization" value="1" scale="1" min="0" max="1000" free="0" />
     <parameter name="Normalization" value="1" scale="1" min="0" max="1000" free="0" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
  </source>

Alternatively the ``DRBPhibarNodes`` model may be used that linearly interpolates for
the background model normalisation for :math:`\bar{\varphi}`-layer that are not present
in the model. Below is an example :ref:`model definition file <glossary_moddef>` for this
model.

.. code-block:: xml

   <source name="Background(0.75-1.0)MeV" type="DRBPhibarNodes" instrument="COM" id="100001">
     <node>
       <parameter name="Phibar"        scale="1.0" value="1.0"  min="0.0" max="50.0"   free="0"/>
       <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
     </node>
     <node>
       <parameter name="Phibar"        scale="1.0" value="3.0"  min="0.0" max="50.0"   free="0"/>
       <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
     </node>
     <node>
       <parameter name="Phibar"        scale="1.0" value="5.0"  min="0.0" max="50.0"   free="0"/>
       <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
     </node>
     <node>
       <parameter name="Phibar"        scale="1.0" value="7.0"  min="0.0" max="50.0"   free="0"/>
       <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
     </node>
     <node>
       <parameter name="Phibar"        scale="1.0" value="9.0"  min="0.0" max="50.0"   free="0"/>
       <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
     </node>
     <node>
       <parameter name="Phibar"        scale="1.0" value="11.0" min="0.0" max="50.0"   free="0"/>
       <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
     </node>
     <node>
       <parameter name="Phibar"        scale="1.0" value="13.0" min="0.0" max="50.0"   free="0"/>
       <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
     </node>
     <node>
       <parameter name="Phibar"        scale="1.0" value="15.0" min="0.0" max="50.0"   free="0"/>
       <parameter name="Normalization" scale="1.0" value="0.0"  min="0.0" max="1000.0" free="0"/>
     </node>
     <node>
       <parameter name="Phibar"        scale="1.0" value="17.0" min="0.0" max="50.0"   free="0"/>
       <parameter name="Normalization" scale="1.0" value="1.0"  min="0.0" max="1000.0" free="1"/>
     </node>
     ...
   </source>

.. warning::
   Depending on the energy band :math:`E'`, a DRB cube may be empty for a number
   of :math:`\bar{\varphi}`-layers, which means that the corresponding
   normalisation cannot be adjusted. To avoid warnings during the :ref:`ctlike`
   model fit it is recommended to fix the ``Normalization`` parameters of the
   corresponding layers. In the first example above, the first four layers are
   fixed. In the second example, all layers with :math:`\bar{\varphi}`
   values below 16 degrees were fixed, the first fitted normalisation
   is for :math:`\bar{\varphi}=17` degrees.
