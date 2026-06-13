# PSBLAS Developer Guide (Internals)

This directory documents the **internal architecture** of PSBLAS: how the
classes relate to each other and how the library is implemented. It is aimed
at developers who modify or extend PSBLAS, not at end users.

The split is deliberate and two-layered:

- **User manual** (`docs/src/*.tex`, built into `psblas-3.9.pdf`) describes the
  public API: what each routine does, its arguments, and the semantics a user
  must know to call it correctly. For example, it documents that `psb_halo`
  accepts a `mode` flag and what synchronous vs. split-phase exchange means for
  the caller.
- **Developer guide** (this directory) describes *how* those features are
  implemented: the dispatch path, the communication-handle class hierarchy, the
  MPI mechanisms behind each scheme, and the extension points. A user never
  needs to read this to use `psb_halo`.

When you add a feature, update the layer that matches its audience. A new public
flag goes in the user manual; a new internal communication scheme goes here.

## Contents

- [communication.md](communication.md) — the communication subsystem: the
  `psb_halo`/`psb_ovrl` → `psi_swapdata` → communication-handle dispatch path,
  the available MPI schemes, and the swap-status state machine.
