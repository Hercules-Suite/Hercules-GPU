#ifndef HERCULES_BUILD_INFO_H
#define HERCULES_BUILD_INFO_H

/* Embedded by Make; independent of the run directory and runtime Git state. */
const char* hercules_build_commit();
const char* hercules_build_tree_state();

#endif
