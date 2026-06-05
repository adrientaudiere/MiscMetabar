# Session Context

## User Prompts

### Prompt 1

claude often ask me (too often) for bash command. How to prevent that behavior

### Prompt 2

# Fewer Permission Prompts

Look through my transcripts' MCP and bash tool calls, and based on those, make a prioritized list of patterns that I should add to my permission allowlist to reduce permission prompts. Focus on read-only commands.

The format for permissions is: `Bash(foo*)`, `Bash(foo)`, `Bash(foo bar *)`, `mcp__slack__slack_read_thread`, etc.

Then, add these to the project `.claude/settings.json` under `permissions.allow`.

## Steps

1. **Locate transcripts.** Session transcript...

