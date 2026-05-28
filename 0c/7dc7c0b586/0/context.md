# Session Context

## User Prompts

### Prompt 1

Base directory for this skill: /home/adrien/.claude/skills/r-pkg-release

Orchestrate an R-package release pipeline. Two modes: **interactive** (default — gated between stages, asks the user before mutating anything) and **`--auto`** (autonomous — runs end-to-end, applies a fixed set of safe deterministic fixes silently, only stops for issues that need human judgment).

## Arguments

`` may contain, in any order:

- A package name (e.g. `MiscMetabar`). If absent, read `Package:` from `DESCRIP...

### Prompt 2

<task-notification>
<task-id>bfmehvbru</task-id>
<tool-use-id>REDACTED</tool-use-id>
<output-file>/tmp/claude-1000/-home-adrien-Nextcloud-IdEst-Projets-pqverse-pqverse-pkg-MiscMetabar/d82706b8-f537-4b1e-98fa-4bc434669314/tasks/bfmehvbru.output</output-file>
<status>completed</status>
<summary>Background command "Run full test suite" completed (exit code 0)</summary>
</task-notification>

### Prompt 3

Base directory for this skill: /home/adrien/.claude/plugins/marketplaces/posit-dev-skills/r-lib/cran-extrachecks

# CRAN Extra Checks

Help R package developers prepare packages for CRAN submission by systematically checking for common ad-hoc requirements that CRAN reviewers enforce but `devtools::check()` doesn't catch.

## Workflow

1. **Initial Assessment**: Ask user if this is first submission or resubmission
2. **Run Standard Checklist**: Work through each item systematically (see below)...

### Prompt 4

Base directory for this skill: /home/adrien/.claude/skills/r-pkg-bump-version

Bump an R package from `X.Y.Z.9xxx` (dev) to a release version `X.Y.Z`. This is a thin wrapper around `usethis::use_version()`. It mutates two files (`DESCRIPTION`, `NEWS.md`) and **does not commit**.

## Arguments

`MiscMetabar patch --auto` may contain, in any order:

- A package name (e.g. `MiscMetabar`). If absent, read `Package:` from `DESCRIPTION` in the current working directory.
- A bump type: `patch` | `mi...

