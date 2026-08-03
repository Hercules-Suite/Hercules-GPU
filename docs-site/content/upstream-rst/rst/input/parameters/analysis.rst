========
Analysis
========
``stiffness_calculation_method``
    * Optional: Yes
    * Type: String
    * Description: The stiffness calculation method. Possible options include ``conventional`` and ``effective``. However, we no longer use ``conventional`` in the recent projects since its performance is concerning, and this parameter is going to be deprecated (the value of it will be set as ``effective`` permanently) for the same reason. This is an optional parameter and the default value is ``effective``.

``use_progressive_meshing``
    * Optional: Yes
    * Type: Boolean
    * Description: Whether to utilize progressive meshing technique. Since it will be more efficient with the utilization of progressive meshing, this parameter will be set as 1 (enabled) permanently and be deprecated in the future. This is an optional parameter and the default value is ``1``.

.. TODO: Actually, from the source code, it seems that this parameter is not a boolean but an integer, and it's the number of the iterations to progressively mesh the model with different frequencies. But I can't see the effect when this number is more than 1 with a minimum working case.

Nonlinear Analysis
==================

(This section is under construction.)

The current version of Hercules can only perform nonlinear analysis with :ref:`Istanbul Velocity Model` enabled.

``include_nonlinear_analysis``
    * Optional: Yes
    * Type: Yes or No
    * Description: Whether to consider nonlinearity in the simulation. This is an optional parameter and the default value is ``no``.

``geostatic_loading_time_sec``
    * Optional: Conditional
    * Type:
    * Description: This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

``geostatic_cushion_time_sec``
    * Optional: Conditional
    * Type:
    * Description: This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

``approximate_geostatic_state``
    * Optional: Conditional
    * Type: Yes or No
    * Description: Whether an initial stress state with :math:`K_o = 1 - \sin\phi` as coefficient of lateral earth pressure at rest should be assumed. This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

``material_plasticity_type``
    * Optional: Conditional
    * Type: String
    * Description: Possible options include ``rate_dependent`` and ``rate_independent``. This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

``material_model``
    * Optional: Conditional
    * Type: String
    * Description: Possible options include ``linear``, ``vonMises_ep``, ``vonMises_fa``, ``vonMises_faM``, ``vonMises_baE``, ``vonMises_baH``, ``vonMises_GQH``, ``vonMises_MKZ``, ``vonMises_RO``, ``MohrCoulomb``, and ``DruckerPrager``. This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

.. Notes: vonMises_FAO and  vonMises_FAM  use kinematic hardening modulus to simulate clays in undrained conditions (very rapid phenomena).
       vonMises_FAO: Original Frederich-Armstrong formulation. Gamma0 is obtained from the modulus reductions curves. Gamma0 defines the maximum linear elastic shear stress “So=G*Gamma0”.
                     psi= Hkin/mu is the ratio between the kinematic hardening modulus and mu.
       vonMises_FAM: Modified Frederich-Armstrong formulation. No elastic region is considered, thus Gamma=0. The modification generates closed hysteric loops under symmetric loading.
       vonMises_EP : Elastoplastic. Gamma0 and psi are set to zero. This model generates elasto-plastic behavior with Tao_max=2Su/sqrt(3).
       vonMises_BAH: Borja and Amies (1994) model for hyperbolic hardening function.  H=psi*mu*(kappa)^m. “psi” and “m” are calibrated from the target G/Gmax curve. “kappa” is a plastic variable
       vonMises_BAE: Borja and Amies (1994) model for exponential hardening function. H= 3mu * kappa^2 / (1+2*kappa). No calibration is needed.
       vonMises_GQH: Groholski et al (2016) Hardening function derivated from ecs. (10) and (13)
       vonMises_MKZ: Matasovic(1994) Hardening function corresponding to the following backbone function: tau_bar=gamma_bar/(1+beta*gamma_bar^s)
       vonMises_RO : Ramberg Osgood. Hardening function corresponding to the following backbone function: gamma=tau/Gmax + alpha * tau^eta/(Gmax*tau_y^(eta-1)).
                     Where phi_RO = tau_y/tao_Max, and tao_Max=2Su/sqrt(3). Variable alpha depends on the offset strain, i.e., alpha = gamma_off*Gmax/tao_Max


``error_tolerance``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

``backbone_errTol``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

``stiff_damp``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

``freq_stiff_damp``
    * Optional: Conditional
    * Type:
    * Description: This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

``tension_cutoff``
    * Optional: Conditional
    * Type: Yes or No
    * Description: Whether to remove tensile stresses or not. This parameter only affects Mohr-Coulomb and Drucker-Prager material models. This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

``material_properties_count``
    * Optional: Conditional
    * Type: Integer
    * Description: This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

``no_substeps``
    * Optional: Conditional
    * Type: Integer
    * Description: This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

``material_properties_list``
    * Optional: Conditional
    * Type: Multiple numbers
    * Format: 
    * Description: This parameter is required when ``include_nonlinear_analysis`` is set to ``yes``.

