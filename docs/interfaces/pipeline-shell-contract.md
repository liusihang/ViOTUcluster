# Pipeline Shell Contract

## Authoritative Boundary

The Python pipeline owns:

- argument validation,
- thread and concurrency budget derivation,
- shell-module directory discovery,
- environment export into shell modules.

Shell modules own:

- external tool invocation,
- per-stage filesystem side effects,
- stage-local logs.

## Required Environment Variables

- `INPUT_DIR`
- `RAW_SEQ_DIR`
- `OUTPUT_DIR`
- `DATABASE`
- `THREADS`
- `THREADS_PER_FILE`
- `MIN_LENGTH`
- `CONCENTRATION_TYPE`
- `Group`
- `ScriptDir`
- `VIOTUCLUSTER_PYTHON`
- `MAX_TASKS`
- `CROSS_VALIDATION_TASKS`
- `TPM_tasks`
- `Assemble_jobs`
- `module_timeout_seconds` is owned by the Python pipeline and enforced before shell modules can hang forever.

## Helper Invocation Rule

Shell modules must invoke Python helpers through package modules:

```bash
"${VIOTUCLUSTER_PYTHON}" -m ViOTUcluster.<module>
```

Shell modules must not assume helper `.py` files are colocated with the shell scripts.
