# Session Context

## User Prompts

### Prompt 1

Base directory for this skill: /home/adrien/.claude/skills/r-pkg-release

Orchestrate an R-package release pipeline. Two modes: **interactive** (default — gated between stages, asks the user before mutating anything) and **`--auto`** (autonomous — runs end-to-end, applies a fixed set of safe deterministic fixes silently, only stops for issues that need human judgment).

## Arguments

`` may contain, in any order:

- A package name (e.g. `MiscMetabar`). If absent, read `Package:` from `DESCRIP...

