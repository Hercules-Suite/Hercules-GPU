====================
Parameter Input File
====================

Basic Rules
===========
When you open a Hercules' parameter input file (default: ``parameters.in``), you will see lines like

.. code-block::

    planes_input_file = inputfiles/parameters.in
    # Istanbul_material_model = yes
    output_mesh = 0
    type_of_damping = bkt2
    region_origin_latitude_deg  = 40.69146389
    number_profile_layers = 8
    profile_layers =
          0     2880   1800    2200     100     100
       1000     3600   2250    2300     100     100
       2000     4000   2500    2350     100     100
       3000     4320   2700    2400     100     100
       4000     4640   2900    2450     100     100
       5000     5770   3490    2600     100     100
      10000     6390   3500    2700     100     100
      20000     6790   3920    2800     100     100

As you can see, the basic rule for writing an parameter input file is

.. code-block::

    parameter_name = parameter_value

or

.. code-block::

    parameter_name =
    parameter_value1 parameter_value2 parameter_value3

A parameter value can be a number (either an integer or floating-point number), a path (of a file or a folder), a string, a boolean value (0 or 1), and a yes/no, etc. A parameter name can also takes multiple values, e.g. ``profile_layers``.

For floating-point number, you can also use scientific notation. E.g., 1.0e-05.

Any line starting with a ``#`` is considered as a comment and will be ignored by Hercules. You can use comments to explain the purpose of the parameter or to disable a parameter without deleting it.

There is no any restriction on how many spaces should be placed between the parameter name, the equal sign, and the parameter value. If you want, you can indent the equal sign or the parameter value based on your preference.

Similarly, there is no restriction on how many spaces should be placed before and among parameter values if the parameter name takes multiple values. However, for parameters that take only one value, spaces are not allowed in the parameter value as they would be considered as delimiters.

Also, you can place parameters in any order in the parameter input file. Just make sure you include all necessary parameters you are going to use.

Parameters
==========

In sections below, we will explain each parameter with following format:

``parameter_name``
    * Optional: Yes/Conditional/No
    * Type: Integer/Float/String/Path/Boolean/Yes or No/Multiple numbers
    * Format: The format for each row below the ``parameter_name``. This column only appears when the type of the parameter is "Multiple numbers".
    * Description:

Hercules would create a new folder for you if the specified path does not exist and the folder is used for outputs. If you notice it does not work for some parameters for outputs, please report it on the `GitHub`_ page.

.. _GitHub: https://github.com/Clifford-Yen/hercules

When Hercules were first developed, it does not accept optional inputs. As a result, most of the parameters explained in sections below are mandatory, and the "Optional" column will only be included for a few of the optional parameters. Another label for the "Optional" column is "Conditional", which means the parameter will be required if you specify a certain value for another parameter. For example, if you specify ``output_displacement = yes``, you will need to specify ``output_displacement_file`` as well. In this case, the "Optional" column for ``output_displacement_file`` will be "Conditional".

In addition, since some parameters are not useful or not really necessary at this point, this documentation will only cover essential parameters. If you can not understand the usage of a parameter and can not find the explanation of it in this documentation, please feel free to leave it as is.

.. toctree::
    :maxdepth: 2
    :caption: You can find parameters used in the parameter input file in the following sections:
 
    material
    dimension
    output
    planewaves
    damping
    source
    time
    topography
    building
    analysis
    drm
