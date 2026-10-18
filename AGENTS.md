# Project Agent Instructions

This repository is a software development project created and maintained under Ben's Agent OS.

Global Agent OS startup is handled by harness-level or system-level instructions. Agents unfamiliar with Agent OS should read `prompts/README.md` before working in this repository.

## Before Working

- Read `prompts/architecture.md` before changing architecture, interfaces, file formats, defaults, or major behavior.
- Read `prompts/plans.md` before continuing implementation, validation, or maintenance work.
- Treat human-authored `prompts/prompt_dev_*.md` files as read-only unless Ben explicitly asks you to edit them.

## Prompt Index

- `prompts/README.md`: Prompt-directory conventions and routing.
- `prompts/architecture.md`: Current architecture, interfaces, formats, and important defaults.
- `prompts/plans.md`: Current state, active work, validation, maintenance guidance, and safe-resume context.

Files matching `prompts/*_deep_records_*.md` are archival. Do not read them during routine development unless Ben explicitly requests historical context or an active context file points to a specific deep record for a current question.

When adding, renaming, or deleting a prompt file, update both this index and `prompts/README.md` in the same change.
