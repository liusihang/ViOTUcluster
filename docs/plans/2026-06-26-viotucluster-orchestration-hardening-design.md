# ViOTUcluster Orchestration Hardening Design

Date: 2026-06-26

## Goal

Harden the Python CLI plus shell-module orchestration layer so that:

1. source-tree and installed layouts share one runtime contract,
2. failed subprocess work exits fast instead of hanging forever,
3. `--disable-binning` preserves the intended viral-only dataflow,
4. thread/concurrency settings reflect the user-visible budget,
5. startup validation fails before expensive downstream work.

## Chosen Approach

Use the existing Python CLI plus shell-module architecture, but tighten the contract instead of rewriting the pipeline:

- keep shell modules as shell modules,
- make Python helpers invokable via `python -m ViOTUcluster.<module>`,
- stop redundant shell-side completion polling after Python workers already block on completion,
- move correctness-sensitive decisions into Python where they are easier to test.

This is lower cost than replacing the orchestration model outright, and it directly addresses the confirmed review issues.

## Planned Changes

### 1. Runtime Contract

- Shell modules will keep `ScriptDir` for locating shell scripts only.
- Python helper entrypoints will be invoked as package modules rather than filesystem-relative `.py` paths.
- The pipeline will export a stable Python executable path for shell modules to reuse.

### 2. Failure Propagation

- Python worker launchers will raise non-zero exit status if any worker future fails.
- Shell wrappers for viral prediction, DRAM, and iPhop will stop performing file-existence polling after the Python driver returns.

### 3. Dataflow Correctness

- `--disable-binning` will stage post-cross-validation viral contigs from `SeprateFile/<sample>/<sample>_filtered.fasta`.

### 4. Resource Budgeting

- Export `THREADS_PER_FILE` from the Python pipeline.
- Replace hard-coded cross-validation parallelism with an environment-driven budget.
- Fix All-in-One default assembly threads to respect `cores // asm_concurrency`.

### 5. Validation

- Call database structure validation during pipeline input checks.

## Test Strategy

- Add focused unit tests for:
  - database structure enforcement,
  - disable-binning source selection,
  - exported orchestration environment variables,
  - assembly-thread default calculation.
- Run standard-library-driven tests where possible because the current local environment does not guarantee `pytest`.
