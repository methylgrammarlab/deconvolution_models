# Plans

Last updated: 2026-10-18

## Current State

- Development is on branch `ben-dev`, created from `main` at commit `367ee0b` (`v0.0.3`).
- The project is a Python package with a Click CLI and model implementations for CelFiE, CelFiE-ISH, ReAtlas, UXM, and Epistate.
- The standard Agent OS project context has been added in `AGENTS.md` and this `prompts/` directory.
- No application code was changed during Agent OS setup.
- The unused Pillow pin has been removed from `requirements.txt`; the codebase contains no Pillow/PIL imports or API calls.

## Active Plan

There is no active implementation request. Before starting feature, bug-fix, refactoring, or dependency work, record the goal, scope, success criteria, and validation commands here.

## Validation Notes

- Confirmed all four required Agent OS files exist and are non-empty.
- Confirmed the root and prompt-directory indexes contain the required active-file and deep-record routing.
- Confirmed the active branch is `ben-dev` and no pre-existing tracked file was modified.
- Checked the new Markdown files for trailing whitespace; none was found.
- Confirmed the source, tests, package metadata, and dependency files contain no remaining Pillow/PIL references and `git diff --check` passes.
- The Python test suite was not run for the unused-dependency cleanup because no runtime code changed.
- Reviewed the four pending project-context Markdown files for publication on 2026-10-18. They contain project documentation only; no application or dependency changes are included.
- Confirmed a single worktree with no nested Git repositories, and that `ben-dev` matched its remote before the documentation commit. No new cache or local-machine files needed ignore rules.
- Checked context-file links, required prompt indexes, and staged whitespace before committing. Runtime tests were not run for this documentation-only change.
- The normal test entry point is `pytest -q`; confirm the required Python environment and `epiread-tools` dependency before relying on that baseline.

## Maintenance Guidance

- Preserve coordinate-system consistency among mixtures, CpG files, intervals, and model-specific atlases.
- Keep CLI behavior, `README.md`, `TUTORIAL.md`, and `prompts/architecture.md` synchronized when options, defaults, formats, or model support change.
- Add focused deterministic tests for numerical behavior. Replace or complete the existing placeholder atlas tests when work touches atlas loading or alignment.
- Keep active context current; archive only strictly historical detail in a dated deep-record file and index it in `prompts/README.md`.

## Safe Resume

1. Confirm the current branch and worktree with `git status --short --branch`.
2. Read `prompts/architecture.md` for the relevant interfaces and constraints.
3. Record the requested work and validation criteria under **Active Plan**.
4. Inspect the affected model wrapper and numerical implementation before editing.
5. Run focused tests, then the broader `pytest -q` suite when the environment supports it.
6. Update this file with completed work, validation results, and remaining risks before handoff.
