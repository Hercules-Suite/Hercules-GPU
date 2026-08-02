========
Overview
========

Hercules is an octree-based finite element earthquake simulator originally developed by the Quake Group at Carnegie Mellon. Hercules is used to solve large-scale earthquake wave propagation problems. For more information and its original source code, please visit `Hercules' GitHub page`_.

.. _Hercules' GitHub page: https://github.com/CMU-Quake/hercules

However, there is nowhere to find usable documentation for Hercules, which makes it really hard to understand and test. This project is aimed to finally document everything on how to use Hercules. The background knowledge and theories used in Hercules are not the priority in this documentation.

Disclaimer
==========
As the maintainer of this project, I, Clifford Yen, am not the owner nor the original developer of Hercules. Our research group has been making improvements to Hercules based on the copy that we got from the original developers and published the source code on `GitHub`_. This documentation is also written based on our version of Hercules. Since it has been evolving for a long time, our version of Hercules is very different from the original version. Please use it at your own risk.

.. _GitHub: https://github.com/Clifford-Yen/hercules

Octree Dimension
================
Since Hercules is an octree-based finite element simulator. All the dimensions in Hercules have to be "octree dimension." That means all dimensions are related to the largest dimension of the domain. Specifically, dimensions have to be 

.. math::
    
    \frac{\text{(The largest dimension of the domain)}}{2^n} \times m

where m and n are positive integers (n is up to 31).

This is why Hercules is so efficient (because all elements are cubical).