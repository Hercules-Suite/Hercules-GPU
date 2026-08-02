==================
Time and Frequency
==================

``simulation_start_time_sec``
    * Optional: Yes
    * Type: Float
    * Description: The start time (in second) for the simulation. This is an optional parameter and the default value is ``0``.

``simulation_end_time_sec``
    * Type: Float
    * Description: The end time (in second) for the simulation.

``simulation_delta_time_sec``
    * Type: Float
    * Description: The time step length (in second) used in the simulation.

``simulation_wave_max_freq_hz``
    * Type: Float
    * Description: The maximum wave frequency (in Hz) for the simulation. 

``simulation_node_per_wavelength``
    * Type: Integer
    * Description: The number of nodes per wavelength.

``simulation_shear_velocity_min``
    * Type: Float
    * Description: The minimum shear velocity. Shear velocities in the simulation that are smaller than this value will be replaced with it.

Since Hercules utilizes explicit FEM, the element size :math:`h` is determined with shear velocity :math:`V_S`, nodes per wavelength :math:`N_p`, and maximum wave frequency :math:`f_{max}`. The formula can be expressed as

.. math::

    h = \frac{V_S}{N_p f_{max}}

And since we also defined the minimum shear velocity :math:`V_{S, min}`, the minimum element size :math:`h_{min}` can be derived by

.. math::

    h_{min} = \frac{V_{S, min}}{N_p f_{max}}

Once you get the minimum element size, the time step :math:`\Delta t` can derive as follows.

.. math::

    \Delta t = \frac{h_{min}}{V_{P, max}}

where :math:`V_{P, max}` is the maximum P wave velocity in the domain.