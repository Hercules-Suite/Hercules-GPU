=============
Miscellaneous
=============

Boolean Handling in Hercules' Source Code
=========================================

For some reason, functions in Hercules' source code that return a boolean or yes/no value do not always use the consistent convention. While some functions return 0 for false and 1 for true, other functions return 0 for true and -1 (or any non-zero value) for false. To make it even more confusing, some functions return 1 as yes and -1 as no.

This inconsistency is not a problem for users who are not intended to modify the Hercules' source code. But for uses who are, be aware of this inconsistency.

Multiple Parameter Input Files
==============================

Hercules used to take more than one ``paramters.in``-like input files. As a result, when digging into Hercules' source code, you may find some comments about "physics.in", "numerical.in", etc. Nevertheless, they are pointing the same and only one parameter input file now.