==========
Topography
==========

``include_topography``
    * Optional: Yes
    * Type: Yes or No
    * Description: Whether to consider topography in the simulation. This is an optional parameter and the default value is ``no``.

``computation_method``
    * Optional: Conditional
    * Type: String
    * Description: Possible options include ``vt`` and ``fem``. This parameter is required if ``include_topography`` is set to ``yes``.

.. TODO: The differences between vt and fem have to be explained.

``maximum_octant_level``
    * Optional: Conditional
    * Type: Integer
    * Description: The maximum octant level allowed. It can be a integer between 1 and 30. As mentioned in :ref:`Time and Frequency`, the minimum element size is calculated as :math:`h_{min} = \frac{V_{S, min}}{N_p f_{max}}`. A proper value for ``maximum_octant_level`` should satisfy the condition :math:`\frac{\text{(The largest dimension of the domain)}}{2^n} < h_{min}`. This parameter is required if ``include_topography`` is set to ``yes``.

.. TODO: ``maximum_octant_level`` should be able to be calculated automatically?

``topographybase_zcoord``
    * Optional: Conditional
    * Type: Float
    * Description: The base z coordinate of the topography. It should be a floating number bigger than ``0`` that satisfy the condition :math:`z_b = \frac{\text{(The largest dimension of the domain)}}{2^n} > e_h`, where :math:`z_b` is ``topographybase_zcoord``, :math:`e_h` is the highest elevation in the domain, and ``n`` is an integer that satisfies the condition. This parameter is required if ``include_topography`` is set to ``yes``. 

``type_of_etree``
    * Optional: Conditional
    * Type: String
    * Description: Possible options include ``sqd`` and ``full``. This parameter is required if ``include_topography`` is set to ``yes``.

.. TODO: The differences between sqd and full have to be explained.

``topography_directory``
    * Optional: Yes
    * Type: Path
    * Description: The path to the directory where the topography information is stored. This is an optional parameter and the default value is ``inputfiles/topography``.

``consider_nonlinear_topography``
    * Optional: Yes
    * Type: Yes or No
    * Description: Whether to consider topography in a nonlinear analysis. This is an optional parameter and the default value is ``no``.