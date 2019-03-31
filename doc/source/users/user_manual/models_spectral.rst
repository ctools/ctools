.. _um_models_spectral:

Spectral model components
-------------------------

The following sections present the spectral model components that are available 
in ctools.

.. warning::

   Source intensities are generally given in units of
   :math:`{\rm photons}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}`.

   **An exception to this rule exists for the** ``DiffuseMapCube`` **spatial
   model where intensities are unitless** and the spectral model presents a
   relative scaling of the diffuse model cube values.

   If spectral models are used in combination with ``RadialAcceptance`` CTA
   background models, intensity units are given in
   :math:`{\rm events}\,\,{\rm s}^{-1}\,{\rm MeV}^{-1}\,{\rm sr}^{-1}`
   and correspond to the on-axis count rate.

   For other CTA background models intensities are unitless and the spectral
   model presents a relative scaling of the background model values.


Constant
^^^^^^^^

  .. code-block:: xml

     <spectrum type="Constant">
       <parameter name="Normalization" scale="1e-16" value="5.7" min="1e-07" max="1000.0" free="1"/>
     </spectrum>

  This spectral model component implements the constant function

  .. math::
     M_{\rm spectral}(E) = N_0

  where

  * :math:`N_0` = ``Normalization``
    :math:`({\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1})`

  .. note::
     For compatibility with the Fermi/LAT ScienceTools the model type
     ``Constant`` can be replaced by ``ConstantValue`` and the parameter
     ``Normalization`` by ``Value``.


Power law
^^^^^^^^^

  .. code-block:: xml

    <spectrum type="PowerLaw">
      <parameter name="Prefactor"   scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
      <parameter name="Index"       scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
      <parameter name="PivotEnergy" scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
    </spectrum>

  This spectral model component implements the power law function

  .. math::
     M_{\rm spectral}(E) = k_0 \left( \frac{E}{E_0} \right)^{\gamma}

  where

  * :math:`k_0` = ``Prefactor``
    :math:`({\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1})`
  * :math:`\gamma` = ``Index``
  * :math:`E_0` = ``PivotEnergy``
    :math:`({\rm MeV})`

  .. warning::
     The ``PivotEnergy`` parameter is not intended to be fitted.

  .. note::
     For compatibility with the Fermi/LAT ScienceTools the parameter
     ``PivotEnergy`` can be replaced by ``Scale``.

  An alternative power law function that uses the integral photon flux as
  parameter rather than the Prefactor is specified by

  .. code-block:: xml

    <spectrum type="PowerLaw">
      <parameter scale="1e-07" name="PhotonFlux" min="1e-07" max="1000.0"    value="1.0" free="1"/>
      <parameter scale="1.0"   name="Index"      min="-5.0"  max="+5.0"      value="-2.0" free="1"/>
      <parameter scale="1.0"   name="LowerLimit" min="10.0"  max="1000000.0" value="100.0" free="0"/>
      <parameter scale="1.0"   name="UpperLimit" min="10.0"  max="1000000.0" value="500000.0" free="0"/>
    </spectrum>

  This spectral model component implements the power law function

  .. math::
     M_{\rm spectral}(E) = \frac{N(\gamma+1)E^{\gamma}}
                                {E_{\rm max}^{\gamma+1} - E_{\rm min}^{\gamma+1}}

  where

  * :math:`N` = ``PhotonFlux``
    :math:`({\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1})`
  * :math:`\gamma` = ``Index``
  * :math:`E_{\rm min}` = ``LowerLimit``
    :math:`({\rm MeV})`
  * :math:`E_{\rm max}` = ``UpperLimit``
    :math:`({\rm MeV})`

  .. warning::
     The ``LowerLimit`` and ``UpperLimit`` parameters are always treated as fixed
     and the flux given by the ``PhotonFlux`` parameter is computed over the
     range set by these two parameters.
     Use of this model allows the errors on the integral flux to be evaluated directly
     by :ref:`ctlike`.

  .. note::
     For compatibility with the Fermi/LAT ScienceTools the model type
     ``PowerLaw`` can be replaced by ``PowerLaw2`` and the parameter
     ``PhotonFlux`` by ``Integral``.


Exponentially cut-off power law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  .. code-block:: xml

    <spectrum type="ExponentialCutoffPowerLaw">
      <parameter name="Prefactor"    scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
      <parameter name="Index"        scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
      <parameter name="CutoffEnergy" scale="1e6"   value="1.0"  min="0.01"  max="1000.0" free="1"/>
      <parameter name="PivotEnergy"  scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
    </spectrum>

  This spectral model component implements the exponentially cut-off power law
  function

  .. math::
     M_{\rm spectral}(E) = k_0 \left( \frac{E}{E_0} \right)^{\gamma}
                           \exp \left( \frac{-E}{E_{\rm cut}} \right)

  where

  * :math:`k_0` = ``Prefactor``
    :math:`({\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1})`
  * :math:`\gamma` = ``Index``
  * :math:`E_0` = ``PivotEnergy``
    :math:`({\rm MeV})`
  * :math:`E_{\rm cut}` = ``CutoffEnergy``
    :math:`({\rm MeV})`

  .. warning::
     The ``PivotEnergy`` parameter is not intended to be fitted.

  .. note::
     For compatibility with the Fermi/LAT ScienceTools the model type
     ``ExponentialCutoffPowerLaw`` can be replaced by ``ExpCutoff`` and
     the parameters ``CutoffEnergy`` by ``Cutoff`` and ``PivotEnergy``
     by ``Scale``.


Super exponentially cut-off power law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  .. code-block:: xml

    <spectrum type="SuperExponentialCutoffPowerLaw">
      <parameter name="Prefactor"    scale="1e-16" value="1.0" min="1e-07" max="1000.0" free="1"/>
      <parameter name="Index1"       scale="-1"    value="2.0" min="0.0"   max="+5.0"   free="1"/>
      <parameter name="CutoffEnergy" scale="1e6"   value="1.0" min="0.01"  max="1000.0" free="1"/>
      <parameter name="Index2"       scale="1.0"   value="1.5" min="0.1"   max="5.0"    free="1"/>
      <parameter name="PivotEnergy"  scale="1e6"   value="1.0" min="0.01"  max="1000.0" free="0"/>
    </spectrum>

  This spectral model component implements the super exponentially cut-off power
  law function

  .. math::
     M_{\rm spectral}(E) = k_0 \left( \frac{E}{E_0} \right)^{\gamma}
                           \exp \left(
                           -\left( \frac{E}{E_{\rm cut}} \right)^{\alpha}
                           \right)

  where

  * :math:`k_0` = ``Prefactor``
    :math:`({\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1})`
  * :math:`\gamma` = ``Index1``
  * :math:`\alpha` = ``Index2``
  * :math:`E_0` = ``PivotEnergy``
    :math:`({\rm MeV})`
  * :math:`E_{\rm cut}` = ``CutoffEnergy``
    :math:`({\rm MeV})`

  .. warning::
     The ``PivotEnergy`` parameter is not intended to be fitted.

  An alternative XML format is supported for compatibility with the Fermi/LAT
  XML format:

  .. code-block:: xml

    <spectrum type="PLSuperExpCutoff">
      <parameter name="Prefactor"   scale="1e-16" value="1.0" min="1e-07" max="1000.0" free="1"/>
      <parameter name="Index1"      scale="-1"    value="2.0" min="0.0"   max="+5.0"   free="1"/>
      <parameter name="Cutoff"      scale="1e6"   value="1.0" min="0.01"  max="1000.0" free="1"/>
      <parameter name="Index2"      scale="1.0"   value="1.5" min="0.1"   max="5.0"    free="1"/>
      <parameter name="Scale"       scale="1e6"   value="1.0" min="0.01"  max="1000.0" free="0"/>
    </spectrum>


Broken power law
^^^^^^^^^^^^^^^^

  .. code-block:: xml

    <spectrum type="BrokenPowerLaw">
      <parameter name="Prefactor"   scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
      <parameter name="Index1"      scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
      <parameter name="BreakEnergy" scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="1"/>
      <parameter name="Index2"      scale="-1"    value="2.70" min="0.01"  max="1000.0" free="1"/>
    </spectrum>

  This spectral model component implements the broken power law function

  .. math::
     M_{\rm spectral}(E) = k_0 \times \left \{
     \begin{eqnarray}
       \left( \frac{E}{E_b} \right)^{\gamma_1} & {\rm if\,\,} E < E_b \\
       \left( \frac{E}{E_b} \right)^{\gamma_2} & {\rm otherwise}
     \end{eqnarray}
     \right .

  where

  * :math:`k_0` = ``Prefactor``
    :math:`({\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1})`
  * :math:`\gamma_1` = ``Index1``
  * :math:`\gamma_2` = ``Index2``
  * :math:`E_b` = ``BreakEnergy``
    :math:`({\rm MeV})`

  .. warning::
     Note that the ``BreakEnergy`` parameter may be poorly constrained if
     there is no clear spectral cut-off in the spectrum.
     This model may lead to complications in the maximum likelihood fitting.

  .. note::
     For compatibility with the Fermi/LAT ScienceTools the parameters
     ``BreakEnergy`` can be replaced by ``BreakValue``.


Smoothly broken power law
^^^^^^^^^^^^^^^^^^^^^^^^^

  .. code-block:: xml

     <spectrum type="SmoothBrokenPowerLaw">
       <parameter name="Prefactor"       scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
       <parameter name="Index1"          scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
       <parameter name="PivotEnergy"     scale="1e6"   value="1.0"  min="0.01"  max="1000.0" free="0"/>
       <parameter name="Index2"          scale="-1"    value="2.70" min="0.01"  max="+5.0"   free="1"/>
       <parameter name="BreakEnergy"     scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="1"/>
       <parameter name="BreakSmoothness" scale="1.0"   value="0.2"  min="0.01"  max="10.0"   free="0"/>
     </spectrum>

  This spectral model component implements the smoothly broken power law function

  .. math::
     M_{\rm spectral}(E) = k_0 \left( \frac{E}{E_0} \right)^{\gamma_1}
                           \left[ 1 +
                           \left( \frac{E}{E_b} \right)^{\frac{\gamma_1 - \gamma_2}{\beta}}
                           \right]^{-\beta}

  where

  * :math:`k_0` = ``Prefactor``
    :math:`({\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1})`
  * :math:`\gamma_1` = ``Index1``
  * :math:`E_0` = ``PivotEnergy``
  * :math:`\gamma_2` = ``Index2``
  * :math:`E_b` = ``BreakEnergy``
    :math:`({\rm MeV})`
  * :math:`\beta` = ``BreakSmoothness``

  .. warning::
     The pivot energy should be set far away from the expected break energy
     value.

  .. warning::
     When the two indices are close together, the :math:`\beta` parameter
     becomes poorly constrained. Since the :math:`\beta` parameter also scales
     the indices, this can cause very large errors in the estimates of the
     various spectral parameters. In this case, consider fixing :math:`\beta`.

  .. note::
     For compatibility with the Fermi/LAT ScienceTools the parameters
     ``PivotEnergy`` can be replaced by ``Scale``,
     ``BreakEnergy`` by ``BreakValue`` and
     ``BreakSmoothness`` by  ``Beta``.


Log parabola
^^^^^^^^^^^^

  .. code-block:: xml

    <spectrum type="LogParabola">
      <parameter name="Prefactor"   scale="1e-17" value="5.878"   min="1e-07" max="1000.0" free="1"/>
      <parameter name="Index"       scale="-1"    value="2.32473" min="0.0"   max="+5.0"   free="1"/>
      <parameter name="Curvature"   scale="-1"    value="0.074"   min="-5.0"  max="+5.0"   free="1"/>
      <parameter name="PivotEnergy" scale="1e6"   value="1.0"     min="0.01"  max="1000.0" free="0"/>
    </spectrum>

  This spectral model component implements the log parabola function

  .. math::
     M_{\rm spectral}(E) = k_0 \left( \frac{E}{E_0} \right)^{\gamma+\eta \ln(E/E_0)}

  where

  * :math:`k_0` = ``Prefactor``
    :math:`({\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1})`
  * :math:`\gamma` = ``Index``
  * :math:`\eta` = ``Curvature``
  * :math:`E_0` = ``PivotEnergy``
    :math:`({\rm MeV})`

  .. warning::
     The ``PivotEnergy`` parameter is not intended to be fitted.

  An alternative XML format is supported for compatibility with the Fermi/LAT
  XML format:

  .. code-block:: xml

     <spectrum type="LogParabola">
       <parameter name="norm"  scale="1e-17" value="5.878"   min="1e-07" max="1000.0" free="1"/>
       <parameter name="alpha" scale="1"     value="2.32473" min="0.0"   max="+5.0"   free="1"/>
       <parameter name="beta"  scale="1"     value="0.074"   min="-5.0"  max="+5.0"   free="1"/>
       <parameter name="Eb"    scale="1e6"   value="1.0"     min="0.01"  max="1000.0" free="0"/>
     </spectrum>

  where

  * ``alpha`` = -``Index``
  * ``beta`` = -``Curvature``


Gaussian
^^^^^^^^

  .. code-block:: xml

     <spectrum type="Gaussian">
       <parameter name="Normalization" scale="1e-10" value="1.0"  min="1e-07" max="1000.0" free="1"/>
       <parameter name="Mean"          scale="1e6"   value="5.0"  min="0.01"  max="100.0"  free="1"/>
       <parameter name="Sigma"         scale="1e6"   value="1.0"  min="0.01"  max="100.0"  free="1"/>
     </spectrum>

  This spectral model component implements the gaussian function

  .. math::
     M_{\rm spectral}(E) = \frac{N_0}{\sqrt{2\pi}\sigma}
                           \exp \left( \frac{-(E-\bar{E})^2}{2 \sigma^2} \right)

  where

  * :math:`N_0` = ``Normalization``
    :math:`({\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1})`
  * :math:`\bar{E}` = ``Mean``
    :math:`({\rm MeV})`
  * :math:`\sigma` = ``Sigma``
    :math:`({\rm MeV})`


File function
^^^^^^^^^^^^^

  .. code-block:: xml

     <spectrum type="FileFunction" file="data/filefunction.txt">
       <parameter scale="1.0" name="Normalization" min="0.0" max="1000.0" value="1.0" free="1"/>
     </spectrum>

  This spectral model component implements an arbitrary function
  that is defined by intensity values at specific energies.
  The energy and intensity values are defined using an ASCII file with
  columns of energy and differential flux values.
  Energies are given in units of
  :math:`{\rm MeV}`,
  intensities are given in units of
  :math:`{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}`.
  The only parameter is a multiplicative normalization:

  .. math::
     M_{\rm spectral}(E) = N_0 \left. \frac{dN}{dE} \right\rvert_{\rm file}

  where

  * :math:`N_0` = ``Normalization``

  .. warning::
     If the file name is given without a path it is expected that the file
     resides in the same directory than the XML file.
     If the file resides in a different directory, an absolute path name should
     be specified.
     Any environment variable present in the path name will be expanded.


Node function
^^^^^^^^^^^^^

  .. code-block:: xml

     <spectrum type="NodeFunction">
       <node>
         <parameter name="Energy"    scale="1.0"   value="1.0" min="0.1"   max="1.0e20" free="0"/>
         <parameter name="Intensity" scale="1e-07" value="1.0" min="1e-07" max="1000.0" free="1"/>
       </node>
       <node>
         <parameter name="Energy"    scale="10.0"  value="1.0" min="0.1"   max="1.0e20" free="0"/>
         <parameter name="Intensity" scale="1e-08" value="1.0" min="1e-07" max="1000.0" free="1"/>
       </node>
     </spectrum>

  This spectral model component implements a generalised broken 
  power law which is defined by a set of energy and intensity values
  (the so called nodes) that are piecewise connected by power laws.
  Energies are given in units of
  :math:`{\rm MeV}`,
  intensities are given in units of
  :math:`{\rm ph}\,\,{\rm cm}^{-2}\,{\rm s}^{-1}\,{\rm MeV}^{-1}`.

  .. warning::
     An arbitrary number of energy-intensity nodes can be defined in a node
     function.
     The nodes need to be sorted by increasing energy.
     Although the fitting of the ``Energy`` parameters is formally possible
     it may lead to numerical complications.
     If ``Energy`` parameters are to be fitted make sure that the ``min``
     and ``max`` attributes are set in a way that avoids inversion of the energy
     ordering.


Composite model
^^^^^^^^^^^^^^^

  .. code-block:: xml

     <spectrum type="Composite">
       <spectrum type="PowerLaw" component="SoftComponent">
         <parameter name="Prefactor"   scale="1e-17" value="3"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"       scale="-1"    value="3.5" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6"   value="1"  min="0.01"  max="1000.0" free="0"/>
       </spectrum>
       <spectrum type="PowerLaw" component="HardComponent">
         <parameter name="Prefactor"   scale="1e-17" value="5"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"       scale="-1"    value="2.0" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6"   value="1"  min="0.01"  max="1000.0" free="0"/>
       </spectrum>
     </spectrum>

  This spectral model component implements a composite model that is the
  sum of an arbitrary number of spectral models, computed using

  .. math::
     M_{\rm spectral}(E) = \sum_{i=0}^{N-1} M_{\rm spectral}^{(i)}(E)

  where :math:`M_{\rm spectral}^{(i)}(E)` is any spectral model component
  (including another composite model), and :math:`N` is the number of
  model components that are combined.


Multiplicative model
^^^^^^^^^^^^^^^^^^^^

  .. code-block:: xml

     <spectrum type="Multiplicative">
       <spectrum type="PowerLaw" component="PowerLawComponent">
         <parameter name="Prefactor"   scale="1e-17" value="1.0"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"       scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6"   value="1.0"  min="0.01"  max="1000.0" free="0"/>
       </spectrum>
       <spectrum type="ExponentialCutoffPowerLaw" component="CutoffComponent">
         <parameter name="Prefactor"    scale="1.0" value="1.0" min="1e-07" max="1000.0" free="0"/>
         <parameter name="Index"        scale="1.0" value="0.0" min="-2.0"  max="+2.0"   free="0"/>
         <parameter name="CutoffEnergy" scale="1e6" value="1.0" min="0.01"  max="1000.0" free="1"/>
         <parameter name="PivotEnergy"  scale="1e6" value="1.0" min="0.01"  max="1000.0" free="0"/>
       </spectrum>
     </spectrum>

  This spectral model component implements a composite model that is the
  product of an arbitrary number of spectral models, computed using

  .. math::
     M_{\rm spectral}(E) = \prod_{i=0}^{N-1} M_{\rm spectral}^{(i)}(E)

  where :math:`M_{\rm spectral}^{(i)}(E)` is any spectral model component
  (including another composite model), and :math:`N` is the number of
  model components that are multiplied.


Exponential model
^^^^^^^^^^^^^^^^^

  .. code-block:: xml

    <spectrum type="Exponential">
      <spectrum type="FileFunction" file="opacity.txt">
	    <parameter scale="-1.0" name="Normalization" min="0.0" max="100.0" value="1.0" free="1"/>
      </spectrum>
    </spectrum>

  This spectral model component implements the exponential of an arbitrary
  spectral model and computes

  .. math::
     M_{\rm spectral}(E) = \exp \left( M_{\rm spectral}(E) \right)

  where :math:`M_{\rm spectral}(E)` is any spectral model component.

  The model can be used to describe a spectrum with EBL absorption based on a
  tabulated model of opacity as a function of photon energy. The corresponding
  XML file structure for such a model is shown below:

  .. code-block:: xml

    <spectrum type="Multiplicative">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor"   scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
        <parameter name="Index"       scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
        <parameter name="PivotEnergy" scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
      </spectrum>
      <spectrum type="Exponential">
        <spectrum type="FileFunction" file="opacity.txt">
          <parameter scale="-1.0" name="Normalization" min="0.0" max="100.0" value="1.0" free="1"/>
        </spectrum>
      </spectrum>
    </spectrum>

  This example corresponds to the function

  .. math::
     M_{\rm spectral}(E) = k_0 \left( \frac{E}{E_0} \right)^{\gamma}
                           \times \exp\left( -\alpha \, \tau(E) \right)

  where
  
  * the first block/factor corresponds to a power law;
  * the second block/factor models EBL absorption, and it points to an
    ASCII file with two columns containing energy in :math:`{\rm MeV}`
    as first column and opacity :math:`\tau`  as second column, respectively;
  * the parameter :math:`\alpha` = ``Normalization`` represents an
    opacity scaling factor.

  .. note::
     The ``Exponential`` model implements the function :math:`y=\exp(x)`,
     hence in the example the ``scale`` attribute of the ``Normalization``
     parameter was set to ``-1`` to implement the form
     :math:`y=\exp(-x)`.
