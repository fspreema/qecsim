# qecsim

**qecsim** is a Python package for simulating quantum error correction codes and lattice surgery operations using [Stim](https://github.com/quantumlib/Stim), a fast stabilizer simulator.

### Features
- Supports multiple quantum error correction codes:
  - Repetition code
  - (Rotated) Surface code
  - XZZX code

- Implements CX gate via lattice surgery on the rotated surface code
  - Following Flows are currently Supported:
    - X -> XX
    - XX -> X
    - x -> X
    - Z -> ZZ
    - ZZ -> Z
    - Z -> Z
    - ZX -> ZX

- Numerical Threshold Approximation
  - Implements numerical threhsold approximations for any given Sinter.TaskStats File

