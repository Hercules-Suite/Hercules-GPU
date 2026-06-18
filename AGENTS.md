# Agent Rules For This Repository

- Do not try to build or run the full Hercules GPU solver locally in this
  workspace. The local machine is not the supported build environment.
- Use local checks only when they do not require the Vista toolchain, such as
  `git diff --check`, code inspection, and targeted text searches.
- Put GPU performance as a priority when implementing everything. Prefer
  approaches that preserve efficient GPU execution, minimize host/device
  transfers, and avoid unnecessary CPU-side work in solver hot paths.
- The user will compile and run the solver on Vista. After code changes, report
  what was changed, what local checks were run, and whether a Vista build/run is
  needed.
- Keep commits focused and push the current branch after committing changes,
  unless the user explicitly says not to push.
