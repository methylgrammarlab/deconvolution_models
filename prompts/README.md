# Project Prompts

This directory holds the active project context that lets humans and agents resume development safely without reconstructing the repository's architecture and implementation state from scratch.

## Standard Files

- `README.md`: Conventions, routing, and the prompt-file index.
- `architecture.md`: Current architecture decisions, project structure, interfaces, file formats, and important defaults.
- `plans.md`: Current implementation state, active plans, validation notes, maintenance guidance, and safe-resume steps.

Optional `prompt_dev_*.md` files contain human-authored implementation instructions. Treat them as read-only unless Ben explicitly asks for edits.

Optional dated `*_deep_records_*.md` files contain historical or archival context. Do not read them during routine development unless Ben explicitly requests historical context or an active file points to a particular deep record for a current question.

## Update Rules

- Keep `architecture.md` and `plans.md` concise, current, and focused on information needed for active development.
- Update `architecture.md` when architecture, interfaces, file formats, defaults, or major behavior changes.
- Update `plans.md` after implementation or maintenance work with what changed, what was validated, and what remains.
- Move strictly historical detail to a dated deep-record file only when archival context is useful; do not archive information still needed for normal maintenance.
- When adding, renaming, or deleting a prompt file, update the prompt indexes in both this file and the root `AGENTS.md`.

## Deep Records Index

Deep-record files are historical archives. Do not read them during routine development unless Ben explicitly asks for historical context or an active file points to one for a current question.

- None yet.
