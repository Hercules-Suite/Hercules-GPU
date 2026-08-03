=============================
3D Velocity Model Input Files
=============================
To make it easier to understand the format of the 3D velocity model input files, here is an example before we convert it to the format that Hercules can read:

.. code-block::
    
    0.0,0.0,0,1709.769,4188.062,1.85
    0.0,0.0,1000,1894.495,4640.546,1.898
    0.0,0.0,3000,2430.156,5952.642,2.02
    0.0,0.0,17400,3400.0,8328.265,2.197
    500.0,0.0,0,1743.876,4271.606,1.859
    500.0,0.0,1000,1920.505,4704.257,1.905
    500.0,0.0,5000,2625.958,6432.257,2.06
    500.0,0.0,17400,3400.0,8328.265,2.197
    0.0,500.0,0,1577.978,3865.241,1.814
    0.0,500.0,1000,1840.914,4509.3,1.885
    0.0,500.0,6000,2714.21,6648.43,2.077
    500.0,500.0,0,1609.918,3943.478,1.823
    500.0,500.0,1000,1861.275,4559.174,1.89
    500.0,500.0,5000,2712.018,6643.06,2.076
    500.0,500.0,17400,3400.0,8328.265,2.197

The meanings of each column are:

    * x (in meter, along the north-south direction with north being positive)
    * y (in meter, along the east-west direction with east being positive)
    * z (in meter, positive downward)
    * Vs (in meter/second)
    * Vp (in meter/second)
    * density (in ton/meter :superscript:`3`)

The order of the columns is not important since this file has to be converted into multiple files which will be explained later. Note that the increments of x and y can be different, but they have to be consistent (i.e. the increment of x has to be the same for all lines, so does the increment of y). The increment of z can also be different, and it can varies from point to point. The number layers can be different for different points, too.

Some rules have to be followed when writing the 3D velocity model input file:

    * The x coordinates should be written in ascending order
    * The y coordinates should be written in ascending order
    * Finish writing all points with the same y coordinate before moving to the next y coordinate.
    * For each x-y coordinate, the z coordinates should be written in ascending order and should be unique.

The file above should be converted into multiple files. The required files include: ``x.in``, ``y.in``, ``z.in``, ``vs.in``, ``vp.in``, ``rho.in``, and ``index.in``. 

The file ``x.in`` contains the x coordinates of all points. The following is the content of ``x.in`` for the example above:

.. code-block::

    2
    0.0
    500.0

The first line of the file is the number of points along the x direction (same as the number of following lines in the file). The following lines are the x coordinates of all points. Note that the x coordinates should be written in ascending order and should be unique. The same rule applies to ``y.in``. In this case, the content of ``y.in`` would be the same as ``x.in``. But if the increment of y or/and the number of points of y are different from x, the content of ``y.in`` will be different.

For ``z.in``, ``vs.in``, ``vp.in``, and ``rho.in``, the content in each of these files is simply the corresponding column in the original file. For example, the content of ``z.in`` for the example above is:

.. code-block::

    15
    0
    1000
    3000
    17400
    0
    1000
    5000
    17400
    0
    1000
    6000
    0
    1000
    5000
    17400

Again, the first line is the total number of the data (same as the number of following lines in the file). The following lines are the z coordinates of all points. The same rule applies to ``vs.in``, ``vp.in``, and ``rho.in``. 

The file ``index.in`` is a special file that includes the starting indices for all points. For example, the content of ``index.in`` for the example above is:

.. code-block::

    4
    0
    4
    8
    11

The first line is the total number of the x-y coordinates (same as the number of the following lines. It's also the same as the multiplication of the number of x points and the number of y points). The following lines are the starting indices for all points. 

For example, the first point (0.0, 0.0) has the starting index of 0 and the second point (500.0, 0.0) has the starting index of 4, which means the first (index 0) to the fourth (index 3) line in ``z.in``, ``vs.in``, ``vp.in``, and ``rho.in`` are the data for the first point. Similarly, the second point (500.0, 0.0) has the starting index of 4 and the third point (0.0, 500.0) has the starting index of 8, which means the fifth (index 4) to the eighth (index 7) line in ``z.in``, ``vs.in``, ``vp.in``, and ``rho.in`` are the data for the second point. So on and so forth.

The conversion can be done easily with the help of Pandas, a Python library. For your convenience, here is a working Python script that can convert the original file to the required files:

.. code-block:: python3

    import pandas as pd

    def convertToHerculesVelocityModel(fileName):
        df = pd.read_csv(fileName)
        for col in ['z', 'vs', 'vp', 'rho']:
            df[col].to_csv(col+'.in', index=False, header=[len(df[col])])
        df = df.drop_duplicates(subset=['x', 'y'])
        df.index.name = len(df)
        indexSeries = df.index.to_series()
        indexSeries.to_csv('index.in', index=False, header=[len(indexSeries)])
        for col in ['x', 'y']:
            df[col].drop_duplicates().to_csv(col+'.in', index=False, header=[len(df[col].drop_duplicates())])

    if __name__ == '__main__':
        convertToHerculesVelocityModel('herculesVelocityModel.csv')

All the converted files should be placed in ``./inputfiles/materialfiles`` in order to let Hercules find them. The original file is not necessary, but it is recommended to keep it to easily understand the meaning of each converted file.
