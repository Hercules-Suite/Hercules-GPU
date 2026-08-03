===========
Plane Waves
===========

.. TEMPLATE
.. ``parameter_name``
..     * Optional: Yes/Conditional/No
..     * Type: Integer/Float/String/Path/Boolean/Yes or No/Multiple numbers
..     * Description:

(This section is under construction.)

Incident Plane Waves
====================

``include_incident_planewaves``
    * Optional: Yes
    * Type: Yes or No
    * Description: This is an optional parameter and the default value is ``no``.

``type_of_wave``
    * Optional: Conditional
    * Type: String
    * Description: Possible options include ``SV`` and ``P``. This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``DRMBox_Noelem_Halfwidth_EW``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``DRMBox_Noelem_Halfwidth_NS``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``DRMBox_Noelem_depth``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``DRMBox_element_size``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``DRM_xo``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``DRM_yo``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``fnc_type``
    * Optional: Conditional
    * Type: String
    * Description: Possible options include ``RICKER`` and ``time_hist``. This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``planewave_strike``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``planewave_Z_angle``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``Ts``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``fc``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``Uo``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``ug_alongstrike``
    * Optional: Conditional
    * Type: Path
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes`` and ``fnc_type = time_hist``.

``ug_alongnormal``
    * Optional: Conditional
    * Type: Path
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``, ``fnc_type = time_hist``, and ``no_components = 2``.

``no_components``
    * Optional: Conditional
    * Type: Integer
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.

``ug_timestep``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_incident_planewaves`` is set to ``yes``.


Homogeneous Half-space
======================

``include_hmg_halfspace``
    * Optional: Yes
    * Type: Yes or No
    * Description: This is an optional parameter and the default value is ``no``.

``etreeBox_xo``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_hmg_halfspace`` is set to ``yes``.

``etreeBox_yo``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_hmg_halfspace`` is set to ``yes``.

``etreeBox_lengthEW``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_hmg_halfspace`` is set to ``yes``.

``etreeBox_lengthNS``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_hmg_halfspace`` is set to ``yes``.

``etreeBox_depth``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_hmg_halfspace`` is set to ``yes``.

``hmgHS_Vs``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_hmg_halfspace`` is set to ``yes``.

``hmgHS_Vp``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_hmg_halfspace`` is set to ``yes``.

``hmgHS_rho``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``include_hmg_halfspace`` is set to ``yes``.
