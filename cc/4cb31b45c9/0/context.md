# Session Context

## User Prompts

### Prompt 1

Base directory for this skill: /home/adrien/.claude/skills/r-audit

Run an external black-box audit of a pqverse R package. The audit creates independent test phyloseq objects, calls every exported function, and produces a re-runnable R script plus a TODO checklist of issues. If a R script is already available here: /home/adrien/Nextcloud/IdEst/Projets/pqverse/audits/{pkg}_script.R propose to pass to phase 4 (run the script).

The target package is specified via `dbpq` (e.g., `/r-audit tidypq...

