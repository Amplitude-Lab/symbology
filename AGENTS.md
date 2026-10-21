# Numerical calculation preferences

The user asks that algorithmic and engineering optimization continue during
calculations, and that measured experience be reused in later calculations.
For high-loop sewing work, consult `docs/hexagon_optimization_lessons.json`
and the current experiment's progress and execution-gate receipts first.

- Profile the actual expensive phase before selecting an optimization. Low
  CPU usage can have different causes; distinguish serial work, memory
  traffic, allocation retention and insufficient task parallelism.
- Preserve completed inputs, coefficients, checkpoints and provenance.
  Validate mathematical changes against independent known results and exact
  reference calculations. Keep assumptions explicit and reserve independent
  physical checks.
- Correctness tests do not establish performance optimality. Describe the
  measured scope and remaining comparisons instead of calling a solver
  fully optimized without evidence.
- The user's latest clarification fixes the >15% CPU qualification threshold
  against V16; comparing a later candidate with V26 does not reset that
  threshold. Retain the 5% RSS tolerance and each measurement's actual version
  identity. Consult the current profile and selection amendments in the
  lessons file; do not confuse benchmark selection with an unverified gate.
- Make physical constraints available to the solver from the start, as
  equations or lazy block generators. Early availability does not require
  early assembly or elimination: schedule blocks by expected benefit, cost,
  fill and memory, and keep expensive blocks pending for later processing.
  Deferral must preserve the equations and affine RHS/origin under coordinate
  changes. Judge scheduling by total CPU-hours within the RSS tolerance,
  not intermediate dimension alone; retain all required final checks.
- Record useful results, unsuccessful alternatives and resource settings in
  the optimization lessons file so future runs do not repeat the same work.
  The current ten-loop job has a 50 GB decimal memory limit; do not infer an
  unlimited memory allowance from low instantaneous usage.
- Package-specific instructions and normalization conventions also apply
  when working inside an extracted calculation package.
