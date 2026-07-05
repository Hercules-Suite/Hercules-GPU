# Repository Guidance

- This is CUDA/MPI GPU code. Keep per-step kernels and time-stepping paths fast;
  prefer setup-time diagnostics over runtime branches, transfers, or extra
  synchronization.
- Preserve CUDA-aware MPI and device-resident synchronization paths. New field
  exchanges should use the existing CUDA-aware helpers when device data is
  already resident.
- Avoid unnecessary local compilation attempts for CUDA changes unless the user
  asks for them; the user compiles and validates on a supercomputer.
- Use focused branches for new work, and commit/push when the user requests it.
- Be explicit in monitor output about physical assumptions, especially where
  acoustic, elastic, and poroelastic coupling paths meet.
