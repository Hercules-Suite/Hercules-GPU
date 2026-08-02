=======
Damping
=======

(This section is under construction.)

.. TODO: Adding a description of different types of damping and their related calculations.

.. NOTE: In damping.h, Hercules treats the damping type as an integer. RAYLEIGH is 0, MASS is 1, NONE is 2, BKT is 3, BKT2 is 4, BKT3 is 5, and BKT3F is 6.
.. As a result, if there is something like ``Param.theTypeOfDamping >= BKT`` in the code, it means that the damping type is BKT, BKT2, BKT3, or BKT3F.

``type_of_damping``
    * Type: String
    * Description: The type of damping Hercules should use. Possible options include ``rayleigh``, ``mass``, ``none``, ``bkt``, ``bkt2``, ``bkt3``, and ``bkt3f``.

``simulation_velocity_profile_freq_hz``
    * Optional: Yes
    * Type: Float
    * Description: This parameter is only effective when ``type_of_damping`` is set to ``bkt``, ``bkt2``, ``bkt3``, or ``bkt3f``. It can be a floating number from ``0.0`` to ``simulation_wave_max_freq_hz``. This is an optional parameter and the default value is ``0.0``.

``use_parametric_q_factor``
    * Optional: Conditional
    * Type: Yes or No
    * Description: This parameter is required if ``type_of_damping`` is set to ``bkt``, ``bkt2``, ``bkt3``, or ``bkt3f``.

``use_infinite_qk``
    * Optional: Conditional
    * Type: Yes or No
    * Description: This parameter is required if ``type_of_damping`` is set to ``bkt``, ``bkt2``, ``bkt3``, or ``bkt3f``.

``parametric_q_factor_constant``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``use_parametric_q_factor`` is set to ``yes``.

``parametric_q_factor_alpha``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``use_parametric_q_factor`` is set to ``yes``.

``parametric_q_factor_beta``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``use_parametric_q_factor`` is set to ``yes``.

``the_threshold_damping``
    * Type: Float
    * Description: The threshold to limit the :math:`\zeta` damping ratio. If the computed :math:`\zeta` is larger than this threshold, it will be set to this threshold.

``do_damping_statistics``
    * Optional: Yes
    * Type: Boolean
    * Description: Whether Hercules should print damping statistics to the standard output (usually the command line window running Hercules). This is an optional parameter and the default value is False (``do_damping_statistics = 0``).

Rayleigh Damping
================
In Hercules, two frequencies used to calculate the Rayleigh damping coefficients are :math:`0.2 f` and :math:`f`, where :math:`f` is the maximum frequency of the simulation (defined by ``simulation_wave_max_freq_hz``). The Rayleigh damping coefficients are calculated as follows:

.. math::
   \frac{\alpha}{\zeta} = \frac{\omega_1 \omega_2 [ (2(\omega_1^2+\omega_2^2) + 2 \omega_1 \omega_2)(\ln(\omega_1) - \ln(\omega_2)) + 3(\omega_2^2 - \omega_1^2)]}{\omega_1^3 - \omega_2^3 + 3 \omega_2^2 \omega_1 - 3 \omega_1^2 \omega_2}

.. math::
   \frac{\beta}{\zeta} = \frac{3 [2 \omega_1 \omega_2 (\ln(\omega_2) - \ln(\omega_1)) + \omega_1^2 - \omega_2^2]}{\omega_1^3 - \omega_2^3 + 3 \omega_2^2 \omega_1 - 3 \omega_1^2 \omega_2}

where :math:`\omega_1 = 2 \pi (0.2 f)` and :math:`\omega_2 = 2 \pi f`. The damping ratio (:math:`\zeta`) is calculated as follows element by element:

.. math::

   \zeta = \frac{25.0}{V_S}

where :math:`V_S` is the shear wave velocity of the element in meters per second. Note that the damping ratio is limited by the threshold damping value (``the_threshold_damping``). Once the damping ratio exceeds this threshold, it will be set to this threshold.


Mass Damping
============
Similar to Rayleigh damping, there are also two frequencies used to calculate the mass damping coefficients: :math:`0.1 f` and :math:`8 f`. The mass damping is only proportional to the mass matrix. As a result, the damping coefficients are calculated as follows:

.. math::
    \frac{\alpha}{\zeta} = 1.3 \frac{2 \omega_1 \omega_2 \ln(\frac{\omega_2}{\omega_1})}{\omega_2 - \omega_1}

.. math::
    \frac{\beta}{\zeta} = 0


.. BKT Damping
.. ===========